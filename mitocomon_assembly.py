#!/usr/bin/env python3

"""
@author: Yoshikazu Furuta

MitoCOMON assembly pipeline

Construct a complete mitochondrial DNA sequence from read output of four amplicon from MinION.

"""


# TODO Write README.md documents.
# TODO Let the whole mitocomon.py to leave the log.

from Bio import SeqIO
from Bio.Seq import Seq
import argparse
import datetime
import os
import subprocess
import sys
import re

version = '3.3.1'

def get_arguments():

	def valid_file(param):  # check if input file ends on .fastq or .fasta
		base, ext = os.path.splitext(param)
		if ext.lower() not in ('.fastq','fq'):
			raise argparse.ArgumentTypeError('File extension must be .fastq or .fq.')
		return param

	def dir_path(string):
		if os.path.exists(string):
			string = os.path.join(os.getcwd(), string)
			return string
		else:
			string = os.path.join(os.getcwd(), string)
			os.makedirs(string) # create the folder
			return string


	parser = argparse.ArgumentParser(description='mitocomon: A pipeline for construction of complete mitochondrial sequences from long PCR fragment sequences' )
	parser.add_argument('-i', '--input', required=True, type = valid_file,
						help='Read file in fastq format. [Required]')
	parser.add_argument('-o', '--outputfolder', required=False, type = dir_path, default='./',
						help='Save the results in the specified outputfolder. Default = current working directory')
	parser.add_argument('-p','--prefix',required=True, type = str,
						help='Prefix of the output file. [Required]')
	parser.add_argument('-a', '--amplicon', required=False, type=str, default="mammal",
						help='Amplicons to analyze. When using for amplicons with custom primers, designate custom and use -c option. Choose from mammal, bird, custom. Default =  mammal')
	parser.add_argument('-c', '--custom_amplicon', required=False, type = valid_file,
						help='File of amplicons with their names and primer pair sequences.')
	parser.add_argument('-t', '--threads', required=False, type=str, default="1",
						help='Number of threads used for each procedure. Default = 1.')
	parser.add_argument('-m','--maxr',required=False, type = str, default=1000,
						help='Max number of reads used in the amplicon_sorter step. Larger number recommended when some fragments showed low concentration, but would take longer calculation time. Default = 1000.')
	parser.add_argument('--medaka_model',required=False, type = str,
						help='The model to use in Medaka. Use if the automatic model detection of Medaka does not work. Check the model list in the help message of medaka_consensus.')
	parser.add_argument('-ldc','--length_diff_consensus',required=False, type=str, default="40",
						help='Length difference consensus parameter for amplicon_sorter. Decrease if you want to avoid nesting of shorter amplicon sequence. Default = 40.')

	args = parser.parse_args()
	return args

def save_arguments(): # save all settings in the parameters.txt file
	outputfolder = args.outputfolder
	with open(os.path.join(outputfolder,'parameters.txt'), 'w') as rf:
		rf.write('-----------------------------------------------------------\n')
		rf.write('mitocomon version: ' + version + '\n')
		rf.write('-----------------------------------------------------------\n')
		rf.write('- date and time = ' + datetime.datetime.now().strftime(
			"%B %d, %Y, %I:%M%p") + '\n')
		rf.write('- input file = ' + args.input + '\n')
		rf.write('- output folder = ' + args.outputfolder + '\n')
		rf.write('- prefix =' + args.prefix + '\n')
		rf.write('- maxr =' + str(args.maxr) + '\n')
		rf.write('- amplicon = ' + args.amplicon + '\n')
		rf.write('- threads = ' + args.threads + '\n')
		if args.medaka_model:
			rf.write('- medaka_model = ' + args.medaka_model + '\n')
		rf.write('-----------------------------------------------------------\n')

if __name__ == '__main__':
	try:
		args = get_arguments()
		infolder_file = args.input
		outputfolder = re.sub('$/','',args.outputfolder)
		infolder, infile = os.path.split(os.path.realpath(infolder_file))
		outprefix = args.prefix
		maxr = args.maxr
		threads = args.threads
		module_dir = os.path.abspath(os.path.dirname(__file__))
		ldc = args.length_diff_consensus
		print(infolder_file,outprefix,args.outputfolder)
		save_arguments() # write all settings in the results.txt file

		print('\n#0: Prepare directories.')
		cmd = 'mkdir -p {0}/01_Readfiltering {0}/02_amplicon_sorter {0}/03_assembly {0}/04_mitos'.format(outputfolder)
		print(cmd)
		res = subprocess.run(cmd,shell=True)


		########################
		# 1 Read filtering    #
		########################
		print('\n#1: Read filtering.')
		print('\n#1-1: Parsing amplicon file.')
		#Get amplicon file path.
		amplicon_file = module_dir+"/amplicon_file/mammal.txt"
		if args.amplicon=="bird":
			amplicon_file = module_dir+"/amplicon_file/bird.txt"
		elif args.amplicon=="custom":
			amplicon_file = args.custom_amplicon

		#Parse amplicon file.
		amplicon_names = []
		amplicon_seqs = []
		amplicon_length = []
		with open(amplicon_file, 'r', encoding='utf-8') as f:
			for line in f:
				line = line.strip()
				columns = line.split('\t')
				amplicon_names.append(columns[0])
				amplicon_seqs.append((columns[1],columns[2]))
				amplicon_length.append(columns[3])

		print('\n#1-2: Separate reads corresponding to each amplicon using Cutadapt.')

		for i in range(len(amplicon_names)):
			linked_adapter1 = amplicon_seqs[i][0]+"..."+Seq(amplicon_seqs[i][1]).reverse_complement()
			linked_adapter2 = amplicon_seqs[i][1]+"..."+Seq(amplicon_seqs[i][0]).reverse_complement()
			cmd = 'cutadapt -O 15 -e 0.2 -j {0} -g ^{1}$ -g ^{2}$ --discard-untrimmed -o {3}/01_Readfiltering/{4}_1.fastq {5} 2>&1 | tee {3}/01_Readfiltering/{4}_cutadapt1_log.txt'.format(threads, linked_adapter1, linked_adapter2, outputfolder, amplicon_names[i], infolder_file)
			print('\nParsing reads of {0}'.format(amplicon_names[i]))
			print(cmd)
			res = subprocess.run(cmd,shell=True)

			#Remove chimeric reads
			cmd = 'cutadapt --action none --discard -O 15 -e 0.2 -j {0} -b {1} -b {2} -o {3}/01_Readfiltering/{4}.fastq {3}/01_Readfiltering/{4}_1.fastq 2>&1 | tee {3}/01_Readfiltering/{4}_cutadapt2_log.txt'.format(threads, amplicon_seqs[i][0], amplicon_seqs[i][1], outputfolder, amplicon_names[i])
			print('\nRemoving chimeric reads.')
			print(cmd)
			res = subprocess.run(cmd,shell=True)

			#Remove tmp fastq file.
			cmd = 'rm {0}/01_Readfiltering/{1}_1.fastq'.format(outputfolder, amplicon_names[i])
			print('\nRemoving tmp fastq file.')
			print(cmd)
			res = subprocess.run(cmd,shell=True)


		print('\n#1-3: Filter each fragment read using chopper.')
		fastqfileprefix = amplicon_names
		for fragment in fastqfileprefix:
			fastqfilepath = outputfolder+"/01_Readfiltering/"+fragment+".fastq"
			if os.path.exists(fastqfilepath):
#				cmd = 'cat {0}/01_Readfiltering/{1}.fastq | chopper -q 10 -l 3000 --headcrop 30 --tailcrop 30 --threads {2} > {0}/01_Readfiltering/{1}_all_q10_l3000.fastq'.format(outputfolder,fragment, threads)
				cmd = 'cat {0}/01_Readfiltering/{1}.fastq | chopper -q 10 -l 3000 --threads {2} > {0}/01_Readfiltering/{1}_all_q10_l3000.fastq'.format(outputfolder,fragment, threads)
				print(cmd)
				res = subprocess.run(cmd,shell=True)

				#Read number check.
				choppedfilepath = outputfolder+"/01_Readfiltering/"+fragment+"_all_q10_l3000.fastq"
				entrylist = list(SeqIO.parse(choppedfilepath,"fastq"))
				num_seq = len(entrylist)

		print('\n#1-4: Prepare a file combining all filtered reads.')
		cmd = 'cat {0}/01_Readfiltering/*_q10_l3000.fastq > {0}/01_Readfiltering/amplicon_reads_all_q10_l3000.fastq'.format(outputfolder)
		print(cmd)
		res = subprocess.run(cmd,shell=True)


		########################
		# 2 amplicon_sorter   #
		########################
		print('\n#2: Amplicon sequence extraction.')
		print('\n#2-1: Getting consensus sequences of each fragment using amplicon_sorter.')
		for i in range(len(fastqfileprefix)):
			fastqfilepath = outputfolder+"/01_Readfiltering/"+fastqfileprefix[i]+"_all_q10_l3000.fastq"
			if(os.path.exists(fastqfilepath)):
#				minl = int(int(amplicon_length[i])*7/10)
				minl = int(int(amplicon_length[i])*8/10)
				maxl = int(int(amplicon_length[i])*13/10)
				cmd = 'python3 {0}/amplicon_sorter.py -i {1}/01_Readfiltering/{2}_all_q10_l3000.fastq -o {1}/02_amplicon_sorter/{2}_as -min {4} -max {5} -np {6} -maxr {3} -ldc {7} 2>&1 | tee {1}/02_amplicon_sorter/{2}_as_log.txt'.format(module_dir,outputfolder,fastqfileprefix[i],maxr,minl,maxl,threads,ldc)
				print(cmd)
				res = subprocess.run(cmd,shell=True)

				#Move temporary files.
				cmd = 'mv {0}/{1}_all_q10_l3000_compare.tmp {0}/{1}_all_q10_l3000_comparelist.pickle amplicon_sorter_2024-10-16.py {0}/02_amplicon_sorter/'.format(outputfolder,fastqfileprefix[i])

		print('\n#2-2: Combine all consensus sequences.')
		cmd = 'cat {0}/02_amplicon_sorter/*_as/*_l3000_consensussequences.fasta > {0}/02_amplicon_sorter/{1}_consensussequences_all.fasta'.format(outputfolder,outprefix)
		print(cmd)
		res = subprocess.run(cmd,shell=True)

		print('\n#2-3: Check results of amplicon_sorter')
		asresultfile = str(outputfolder+"/02_amplicon_sorter/"+outprefix+"_consensussequences_all.fasta")
		amplicons = list(SeqIO.parse(asresultfile,'fasta'))
		num_amplicons = len(amplicons)
		if (num_amplicons<len(fastqfileprefix)):
			print("\nOnly {0} fragment sequences were found, less than expected {1}.".format(num_amplicons, len(fastqfileprefix)))
			print("\nContinue anyway.")
		with open(os.path.join(outputfolder,'parameters.txt'), 'a') as rf:
			rf.write("Fragments detected by amplicon_sorter: {0}".format(num_amplicons))


		###################################
		# 3 Minimap2 + Miniasm + Medaka   #
		###################################
		print('\n#3: Assembly and correction.')
		print('\n#3-1: Assemble amplicon sequences using minimap2 and miniasm.')
		cmd = 'minimap2 -x ava-ont {0}/02_amplicon_sorter/{1}_consensussequences_all.fasta {0}/02_amplicon_sorter/{1}_consensussequences_all.fasta > {0}/03_assembly/{1}_consensussequences_all.paf'.format(outputfolder,outprefix)
		print(cmd)
		res = subprocess.run(cmd,shell=True)

		#Filter paf file.
		paffilepath = outputfolder+"/03_assembly/"+outprefix+"_consensussequences_all.paf"
		pafoutfilepath = outputfolder+"/03_assembly/"+outprefix+"_consensussequences_all_filtered.paf"
		seqdivpattern = re.compile(r"dv:f:(.*?)\t")
		idpattern = re.compile(r"consensus_(.*?)_all_q10_l3000")
		with open(paffilepath,'r') as f:
			with open(pafoutfilepath,'w') as of:
				for line in f:
	#				line = line.strip()
					match = seqdivpattern.search(line)
					if match and float(match.group(1))<0.01:
#					if match:
						col = line.split('\t')
						matchid1 = idpattern.search(col[0])
						matchid2 = idpattern.search(col[5])
						if str(matchid1.group(1))!=str(matchid2.group(1)) and \
							(int(col[2])<100 or (int(col[1])-int(col[3]))<100) and \
							(int(col[7])<100 or (int(col[6])-int(col[8]))<100):
							of.write(line)

		cmd = 'miniasm -e 0 -s 300 -o 300 -n 0 -1 -2 -f {0}/02_amplicon_sorter/{1}_consensussequences_all.fasta {0}/03_assembly/{1}_consensussequences_all_filtered.paf > {0}/03_assembly/{1}_consensussequences_assembled.gfa'.format(outputfolder,outprefix)
		print(cmd)
		res = subprocess.run(cmd,shell=True)

		cmd = "awk '/^S/{{print \">\" $2 \"\\n\" $3}}' {0}/03_assembly/{1}_consensussequences_assembled.gfa > {0}/03_assembly/{1}_consensussequences_all.circularise.fasta".format(outputfolder,outprefix)
		print(cmd)
		res = subprocess.run(cmd,shell=True)

		print('\n#3-2: Sequence correction using Medaka.')
		if args.medaka_model:
			cmd = 'medaka_consensus -i {0}/01_Readfiltering/amplicon_reads_all_q10_l3000.fastq -d {0}/03_assembly/{1}_consensussequences_all.circularise.fasta -o {0}/03_assembly/medaka -t {3} -m {2} && \
				medaka_consensus -i {0}/01_Readfiltering/amplicon_reads_all_q10_l3000.fastq -d {0}/03_assembly/medaka/consensus.fasta -o {0}/03_assembly/medaka2 -t {3} -m {2} && \
				medaka_consensus -i {0}/01_Readfiltering/amplicon_reads_all_q10_l3000.fastq -d {0}/03_assembly/medaka2/consensus.fasta -o {0}/03_assembly/medaka3 -t {3} -m {2} && \
				mv {0}/03_assembly/medaka3/consensus.fasta {0}/03_assembly/{1}_consensussequences_all.circularise.medaka.fasta'.format(outputfolder,outprefix,args.medaka_model,threads)
		else:
			cmd = 'medaka_consensus -i {0}/01_Readfiltering/amplicon_reads_all_q10_l3000.fastq -d {0}/03_assembly/{1}_consensussequences_all.circularise.fasta -o {0}/03_assembly/medaka -t {2} && \
				medaka_consensus -i {0}/01_Readfiltering/amplicon_reads_all_q10_l3000.fastq -d {0}/03_assembly/medaka/consensus.fasta -o {0}/03_assembly/medaka2 -t {2} && \
				medaka_consensus -i {0}/01_Readfiltering/amplicon_reads_all_q10_l3000.fastq -d {0}/03_assembly/medaka2/consensus.fasta -o {0}/03_assembly/medaka3 -t {2} && \
				mv {0}/03_assembly/medaka/consensus.fasta {0}/03_assembly/{1}_consensussequences_all.circularise.medaka.fasta'.format(outputfolder,outprefix,threads)
		print(cmd)
		res = subprocess.run(cmd,shell=True,executable='/bin/bash')

		print('\n#3-3: Rotate the start position according to the direction and position of tRNA-Phe.')
		cmd = 'python3 {0}/setstartposition_mtDNA.py {1}/03_assembly/{2}_consensussequences_all.circularise.medaka.fasta {1}/{2}_mtDNAcontig'.format(module_dir,outputfolder,outprefix)
		print(cmd)
		res = subprocess.run(cmd,shell=True)


		############################
		# 4 Annotation by MITOS2   #
		############################
		print('\n#4: Annotation using MITOS2.')
		print('\n#4-1: Preparing reference')
		cmd = 'wget -O {0}/04_mitos/refseq89m.tar.bz2 https://zenodo.org/records/4284483/files/refseq89m.tar.bz2?download=1 && \
			tar -jxf {0}/04_mitos/refseq89m.tar.bz2 -C {0}/04_mitos'.format(outputfolder)
		print(cmd)
		res = subprocess.run(cmd,shell=True)

		print('\n#4-2: Running MITOS2')
		contigfilepath = outputfolder+"/"+outprefix+"_mtDNAcontig.fasta"
		trnalist = ["trnF(gaa)","trnV(tac)","trnL2(taa)","trnI(gat)","trnQ(ttg)","trnM(cat)","trnW(tca)","trnA(tgc)","trnN(gtt)","trnC(gca)","trnY(gta)","trnS2(tga)","trnD(gtc)","trnK(ttt)","trnG(tcc)","trnR(tcg)","trnH(gtg)","trnS1(gct)","trnL1(tag)","trnE(ttc)","trnT(tgt)","trnP(tgg)"]
		rrnalist = ["rrnS","rrnL"]
		proteinlist = ["nad1","nad2","cox1","cox2","atp8","atp6","cox3","nad3","nad4l","nad4","nad5","nad6","cob"]

		for record in SeqIO.parse(contigfilepath, "fasta"):
			cmd = 'mkdir -p {0}/04_mitos/{1}'.format(outputfolder,record.id)
			print(cmd)
			res = subprocess.run(cmd,shell=True)

			out_file = outputfolder+"/04_mitos/"+record.id+"/"+outprefix+"_mtDNA"+record.id+".fasta"
			SeqIO.write(record, out_file, "fasta")

			cmd = 'runmitos.py -i {0}/04_mitos/{1}/{2}_mtDNA{1}.fasta -c 2 -R {0}/04_mitos -r refseq89m -o {0}/04_mitos/{1}/ --noplots'.format(outputfolder,record.id,outprefix)
			print(cmd)
			res = subprocess.run(cmd,shell=True)

			#Change filenames.
			cmd = 'mv {0}/04_mitos/{1}/result.bed {0}/04_mitos/{1}/{2}_mtDNA{1}.bed && \
					mv {0}/04_mitos/{1}/result.faa {0}/04_mitos/{1}/{2}_mtDNA{1}.faa && \
					mv {0}/04_mitos/{1}/result.fas {0}/04_mitos/{1}/{2}_mtDNA{1}.fas && \
					mv {0}/04_mitos/{1}/result.geneorder {0}/04_mitos/{1}/{2}_mtDNA{1}.geneorder && \
					mv {0}/04_mitos/{1}/result.gff {0}/04_mitos/{1}/{2}_mtDNA{1}.gff && \
					mv {0}/04_mitos/{1}/result.mitos {0}/04_mitos/{1}/{2}_mtDNA{1}.mitos && \
					mv {0}/04_mitos/{1}/result.png {0}/04_mitos/{1}/{2}_mtDNA{1}.png && \
					mv {0}/04_mitos/{1}/result.seq {0}/04_mitos/{1}/{2}_mtDNA{1}.seq\
					'.format(outputfolder,record.id,outprefix)
#			print(cmd)
			res = subprocess.run(cmd,shell=True)

			#Completeness check
			bedfilepath = outputfolder+"/04_mitos/"+record.id+"/"+outprefix+"_mtDNA"+record.id+".bed"
			genesfound = []
			with open(bedfilepath,"r") as f:
				for line in f:
					line = line.strip()
					columns = line.split('\t')
					genesfound.append(columns[3])
			#tRNA check
			missingtrna = []
			for trna in trnalist:
				trna_parts = trna.split('(')
				pattern = re.compile(trna_parts[0])
				matchflag = 0
				for gene in genesfound:
					if pattern.match(gene):
						matchflag = 1
				if matchflag==0:
					missingtrna.append(trna)
			#rRNA check
			missingrrna = []
			for rrna in rrnalist:
				pattern = re.compile(rrna)
				matchflag = 0
				for gene in genesfound:
					if pattern.match(gene):
						matchflag = 1
				if matchflag == 0:
					missingrrna.append(rrna)
			#protein check
			missingprotein = []
			for protein in proteinlist:
				pattern = re.compile(protein)
				matchflag = 0
				for gene in genesfound:
					if pattern.match(gene):
						matchflag = 1
				if matchflag == 0:
					missingprotein.append(protein)

			checkfilepath = outputfolder+"/04_mitos/"+record.id+"/genecheck.txt"
			with open(checkfilepath,"w") as af:
				af.write("Check of genes detected in "+record.id+"\n")
				af.write("tRNA missing:")
				for trna in missingtrna:
					af.write(trna)
				af.write("\nrRNA missing:")
				for rrna in missingrrna:
					af.write(rrna)
				af.write("\nProtein missing:")
				for protein in missingprotein:
					af.write(protein)

		print('\nALL DONE!')

	except KeyboardInterrupt:
		sys.exit()