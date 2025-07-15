#!/usr/bin/env python3

"""
@author: Yoshikazu Furuta

MitoCOMON design pipeline

Extract primer candidate sequences which is specific to the taxonomy of interest.

"""


from Bio import Seq, SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqUtils import IUPACData
from collections import Counter
import math
import argparse
import datetime
import os
import subprocess
import sys
import re
import random
from subprocess import PIPE
import statistics
import glob

version = '1.0.0'

#----------------------------------------------------------------------------------------------------
def get_arguments():
	def valid_fastafile(param):  # check if input file ends on .fastq or .fasta
		base, ext = os.path.splitext(param)
		if ext.lower() not in ('.fasta','.fa','.fna'):
			raise argparse.ArgumentTypeError('File extension must be .fasta, .fa, or .fna.')
		return param

	def dir_path(string):
		if os.path.exists(string):
			string = os.path.join(os.getcwd(), string)
			return string
		else:
			string = os.path.join(os.getcwd(), string)
			os.makedirs(string) # create the folder
			return string

	parser = argparse.ArgumentParser(description='Designing primer sets for mitocommon.' )
	parser.add_argument('-i', '--input', required=True, type = valid_fastafile,
						help='Fasta file of Mitochondria sequences. Recommended to download from RefSeq release. [Required]')
	parser.add_argument('-o', '--outputfolder', required=False, type = dir_path, default='./',
						help='Save the results in the specified outputfolder. Default = current working directory')
	parser.add_argument('-p','--prefix',required=True, type = str,
						help='Prefix of the output file. [Required]')
	parser.add_argument('-t','--taxid',required=True, type = str,
						help='Taxon ID of the organism of your interest (e.g. Mammalia, Aves, etc.). You can search it at NCBI Taxonomy or use taxonkit name2taxid. [Required]')
	parser.add_argument('-d', '--datadir',required=True, type= str,
						help='Path of nucl_gb.accession2taxid. Download from https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz [Required]')
	parser.add_argument('-z', '--taxkitoption',required=False, type= str, default = "",
						help='Additional options for running taxonkit. For example, if you keep taxdump.tar.gz in the specific folder, you can add "--data-dir ./PATH/TO/TAXDUMP.tar.gz". Optional.')
	parser.add_argument('-s', '--entropythreshold',required=False, type= float, default = 1.80,
						help='Threshold of shannon entropy score for determining conserved bases. Default is 1.80. Optional.')
	parser.add_argument('-b', '--baseposRefSeqid',required=False, type= str, default = "",
						help='RefSeq ID of the entry you want to use for showing the corresponding position of primer sequences. If not set, the first RefSeq ID in 01_sequence_splitting/tmp/target_refseqid.txt will be used. Optional.')
	parser.add_argument('-e', '--ppenv',required=True, type= str,
						help='The name of Python 2 environment where PrimerProspector is installed. [Required]')

	args = parser.parse_args()
	return args

#----------------------------------------------------------------------------------------------------

def setstartposition(outputfolder, fastafile_base):
	startpos = 0
	endpos = 0
	strand = 0
	score = 0
	infile = outputfolder+"/01_sequence_splitting/"+fastafile_base+".fasta"
	outfile = fastafile_base+"_setstart.fasta"

	with open(outfile, 'w', encoding='utf-8') as f:
		with open(infile,"r") as fasta_file:
			for seq_record in SeqIO.parse(fasta_file,"fasta"):
				with open('./tmptrnascantarget.fasta','w') as ftmp:
					ftmp.write(seq_record.format('fasta-2line'))
				path = os.path.abspath('./tmptrnascantarget.fasta')
				cmd = f'tRNAscan-SE {path} -M vert -q'
				res = subprocess.Popen(cmd,shell=True,stdout=PIPE)

				#Store tRNAscan-SE results in a list.
				results = str(res.stdout.read()).split('\\n')

				score = 0
				strand = 0

				for result in results:
					factors = result.split('\\t')
					if (len(factors)<9):
						continue
					if factors[4] == 'Phe' and float(factors[8])>float(score):
						startpos = int(factors[2])
						endpos = int(factors[3])
						if endpos-startpos > 0:
							strand = 1
						else:
							strand = -1

				#If tRNA-Phe was not found, rotate the sequence and retry.
				if strand == 0:
					rotatedseq = seq_record.seq[999:len(seq_record.seq)] + seq_record.seq[0:999]
					with open('./tmptrnascantarget2.fasta','w') as ftmp2:
						ftmp2.write(f'>{seq_record.description} 1k_rotated\n{seq_record.seq}')
					path = os.path.abspath('./tmptrnascantarget2.fasta')
					cmd = f'tRNAscan-SE {path} -M vert -q'
					res = subprocess.Popen(cmd,shell=True,stdout=PIPE)

					#Store tRNAscan-SE results in a list.
					results = str(res.stdout.read()).split('\\n')

					score = 0
					strand = 0

					for result in results:
						factors = result.split('\\t')
						if (len(factors)<9):
							continue
						if factors[4] == 'Phe' and float(factors[8])>float(score):
							startpos = int(factors[2])
							endpos = int(factors[3])
							if endpos-startpos > 0:
								strand = 1
							else:
								strand = -1
					startpos = startpos + 999
					if (startpos>len(seq_record.seq)):
						startpos = startpos - len(seq_record.seq)

				rotated_seq = ''
				if strand == 1:
					subseq1 = seq_record.seq[0:startpos-1]
					subseq2 = seq_record.seq[startpos-1:len(seq_record.seq)]
					rotated_seq = subseq2+subseq1
					f.write(f'>{seq_record.description} rotated_{startpos}\n{rotated_seq}\n')
				elif strand == -1:
					subseq1 = seq_record.seq[0:endpos]
					subseq2 = seq_record.seq[endpos:len(seq_record.seq)]
					rotated_seq = subseq1.reverse_complement()+subseq2.reverse_complement()
					f.write(f'>{seq_record.description} rotated_c{startpos}\n{rotated_seq}\n')
	cmd = 'rm tmptrnascantarget.fasta tmptrnascantarget2.fasta'
	res = subprocess.run(cmd,shell=True)

#----------------------------------------------------------------------------------------------------

def entropy_calc(target_aligned_file, outputname, outfilepath, threshold):
	align = AlignIO.read(target_aligned_file,"fasta")
	alignmentlength = align.get_alignment_length()
	samplenum = 0
	output_ind = 0
	for i in range(len(align)):
		if (re.match(r'Consensus',align[i].id)==None):
			samplenum = samplenum + 1
		if (re.match(outputname,align[i].id)):
			output_ind = i
	poscount = 0
	poslist = []
	alignmentpos = []
	entlist = []
	for i in range(alignmentlength):
		col = align[:,i]
		if(col[output_ind]=='a' or col[output_ind]=='t' or col[output_ind]=='g' or col[output_ind]=='c'):
			poscount = poscount + 1
			dk = {'a':0,'t':0,'g':0,'c':0,'o':0}
			for j in range(samplenum):
				if(col[j]=='a' or col[j]=='t' or col[j]=='g' or col[j]=='c'):
					dk[col[j]] = dk[col[j]]+1
				else:
					dk['o'] = dk['o'] + 1
			total = dk['a']+dk['t']+dk['g']+dk['c']
			ent = 0;
			if (dk['a']!=0):
				ent = ent - (dk['a']/total)*(math.log2(dk['a']/total))
			if (dk['t']!=0):
				ent = ent - (dk['t']/total)*(math.log2(dk['t']/total))
			if (dk['g']!=0):
				ent = ent - (dk['g']/total)*(math.log2(dk['g']/total))
			if (dk['c']!=0):
				ent = ent - (dk['c']/total)*(math.log2(dk['c']/total))
			if (dk['o']>=0.1*samplenum):
				ent = 2;
			poslist.append(poscount)
			alignmentpos.append(i)
			entlist.append(2-ent)
	#Print the entropy list to a file.
	entropyinfo_list = []
	with open(outfilepath,'w',encoding='utf-8') as f:
		tmppos = 0
		while tmppos < poscount-20:
			tt = 0
			for j in range(tmppos,tmppos+20):
				tt = tt + entlist[j]
			tt = tt/20
			f.write(align[output_ind].id+"\t"+str(tmppos+1)+"\t"+str(tmppos+20)+"\t"+str(alignmentpos[tmppos+1])+"\t"+str(alignmentpos[tmppos+20])+"\t"+str(tt)+"\n")
			if (tt>=threshold):
				entropyinfo_list.append((int(tmppos+1),int(tmppos+20),float(tt),align[output_ind].id,int(alignmentpos[tmppos+1]),int(alignmentpos[tmppos+20])))
			tmppos = tmppos + 1
	return entropyinfo_list

#----------------------------------------------------------------------------------------------------

def extractconservedregionseq(entropyinfo_list,outfilename):
	positions_list = entropyinfo_list
	prevend = 0
	maxstart = 0
	maxend = 0
	maxscore = 0
	rstart = 0
	alprevend = 0
	almaxstart = 0
	almaxend = 0
	alrstart = 0
	return_list = []
	with open(outfilename,'w',encoding='utf-8') as f:
		f.write("specname"+"\t"+"region_start"+"\t"+"region_end"+"\t"+"alignment_start"+"\t"+"alignment_end"+"\t"+"best_region_start"+"\t"+"best_region_end"+"\t"+"best_alignment_start"+"\t"+"best_alignment_end"+"\t"+"information_content"+"\n")
		for start,end,score,specname,alstart,alend in positions_list:
			if prevend==0:
				prevend = end
				alprevend = alend
				maxstart = start
				almaxstart = alstart
				maxend = end
				almaxend = alend
				maxscore = score
				rstart = start
				alrstart = alstart
				continue
			if start > prevend:
				maxscore = round(float(maxscore),4)
				f.write(specname+"\t"+str(rstart)+"\t"+str(prevend)+"\t"+str(alrstart)+"\t"+str(alprevend)+"\t"+str(maxstart)+"\t"+str(maxend)+"\t"+str(almaxstart)+"\t"+str(almaxend)+"\t"+str(maxscore)+"\n")
				return_list.append((almaxstart, almaxend, maxstart, maxend, specname))
				maxstart = start
				almaxstart = alstart
				maxend = end
				almaxend = alend
				maxscore = score
				rstart = start
				alrstart = alstart
			if score > maxscore:
				maxstart = start
				almaxstart = alstart
				maxend = end
				almaxend = alend
				maxscore = score
			prevend = end
			alprevend = alend
		maxscore = round(float(maxscore),4)
		f.write(specname+"\t"+str(rstart)+"\t"+str(prevend)+"\t"+str(alrstart)+"\t"+str(alprevend)+"\t"+str(maxstart)+"\t"+str(maxend)+"\t"+str(almaxstart)+"\t"+str(almaxend)+"\t"+str(maxscore)+"\n")
		return_list.append((almaxstart, almaxend, maxstart, maxend, specname))
	return return_list

#----------------------------------------------------------------------------------------------------

def getconsensus_single(start,end,alignment):
	seqnum = len(alignment)

	# Initialize a dictionary to count bases at each position
	base_counts = {i: Counter() for i in range(start, end)}

	# Iterate through each sequence in the alignment
	for record in alignment:
		seq = record.seq[start:end]
		for i, base in enumerate(seq):
			base_counts[start + i][base] += 1

	# Calculate the consensus sequence with degenerate bases allowed
	consensus_sequence = Seq.Seq("")

	for position in range(start, end):
		consbase = "N"
		counts = base_counts[position]
		most_common_base, most_common_count = counts.most_common()[0]
		if most_common_base=="-" and float(most_common_count/seqnum) > 0.5:
			consbase = "-"
		else:
			if float(most_common_count/seqnum) > 0.8:
				consbase = most_common_base
			elif len(counts.most_common()) > 1:
				secondbest_base, secondbest_count = counts.most_common()[1]
				if secondbest_base == "-":
					if len(counts.most_common()) > 2:
						secondbest_base = counts.most_common()[2][0]
						secondbest_count = counts.most_common()[2][1]
					else:
						consbase = most_common_base
						consensus_sequence += consbase
						continue
				if float((most_common_count+secondbest_count)/seqnum) > 0.8:
					if most_common_base.upper() == "A" and secondbest_base.upper() == "T":
						consbase = "W"
					elif most_common_base.upper() == "T" and secondbest_base.upper() == "A":
						consbase = "W"
					elif most_common_base.upper() == "G" and secondbest_base.upper() == "C":
						consbase = "S"
					elif most_common_base.upper() == "C" and secondbest_base.upper() == "G":
						consbase = "S"
					elif most_common_base.upper() == "A" and secondbest_base.upper() == "G":
						consbase = "R"
					elif most_common_base.upper() == "G" and secondbest_base.upper() == "A":
						consbase = "R"
					elif most_common_base.upper() == "T" and secondbest_base.upper() == "C":
						consbase = "Y"
					elif most_common_base.upper() == "C" and secondbest_base.upper() == "T":
						consbase = "Y"
					elif most_common_base.upper() == "A" and secondbest_base.upper() == "C":
						consbase = "M"
					elif most_common_base.upper() == "C" and secondbest_base.upper() == "A":
						consbase = "M"
					elif most_common_base.upper() == "G" and secondbest_base.upper() == "T":
						consbase = "K"
					elif most_common_base.upper() == "T" and secondbest_base.upper() == "G":
						consbase = "K"
		if(consbase!="-"):
			consensus_sequence += consbase.upper()
	return consensus_sequence

#----------------------------------------------------------------------------------------------------

def getconsensus_fromalignment(regionlist, alignment_file):
	return_list = []

	# Define the input multiple alignment file and the positions to extract
	positions_list = regionlist

	# Read the multiple alignment file
	alignment = AlignIO.read(alignment_file, "fasta")

	for start_position,end_position,start,end,specname in positions_list:
		start_position = int(start_position)
		end_position = int(end_position)
		#Modify the start_position as sequence position of biopython is 0-based.
		start_position = start_position-1
		consensus_sequence = getconsensus_single(start_position,end_position,alignment)
		consensus_sequence_up10 = getconsensus_single(start_position-10,start_position,alignment)
		consensus_sequence_down10 = getconsensus_single(end_position,end_position+10,alignment)

		return_list.append((specname,start,end,consensus_sequence,consensus_sequence_up10,consensus_sequence_down10,start_position+1,end_position))
#		print(specname,start,end,consensus_sequence,consensus_sequence_up10,consensus_sequence_down10,start_position+1,end_position,sep='\t')
	return return_list

#----------------------------------------------------------------------------------------------------

def expand_degenerate_bases(sequence):
    def expand_recursive(seq, pos):
        if pos == len(seq):
            return [seq]
        base = seq[pos]
        if base in IUPACData.ambiguous_dna_values:
            expanded_seqs = []
            for possible_base in IUPACData.ambiguous_dna_values[base]:
                new_seq = seq[:pos] + possible_base + seq[pos+1:]
                expanded_seqs.extend(expand_recursive(new_seq, pos+1))
            return expanded_seqs
        else:
            return expand_recursive(seq, pos+1)
    return expand_recursive(sequence.upper(), 0)

#----------------------------------------------------------------------------------------------------

def run_primer3_seqlist(seqid,sourceseq,seqlist,primer3_inputfile_name,primer3_outputfile_name):
	#Make inputfile for primer3.
	with open(primer3_inputfile_name,"w") as file:
		for i, seq in enumerate(seqlist):
			file.write(f"SEQUENCE_ID={seqid}_{i}\n")
			file.write(f"SEQUENCE_PRIMER={seq}\n")
			file.write(f"PRIMER_TASK=check_primers\n")
			file.write(f"PRIMER_MIN_TM=56\n")
			file.write(f"PRIMER_MAN_TM=64\n")
			file.write(f"PRIMER_MIN_GC=40\n")
			file.write(f"PRIMER_MAN_GC=60\n")
			file.write(f"PRIMER_MAX_POLY_X=5\n")
			file.write(f"PRIMER_MAX_END_STABILITY=4\n")
			file.write(f"PRIMER_TM_FORMULA=1\n")
			file.write(f"PRIMER_SALT_CORRECTIONS=1\n")
			file.write(f"PRIMER_SALT_MONOVALENT=50.0\n")
			file.write(f"PRIMER_SALT_DIVALENT=2.5\n")
			file.write(f"PRIMER_DNTP_CONC=0.3\n")
			file.write(f"PRIMER_MAX_SELF_ANY_TH=47.00\n")
			file.write(f"PRIMER_MAX_SELF_END_TH=47.00\n")
			file.write(f"PRIMER_MAX_HAIRPIN_TH=47.00\n")
			file.write(f"PRIMER_PICK_ANYWAY=1\n")
			file.write(f"=\n")

	#Run Primer3.
	command = f"primer3_core < {primer3_inputfile_name} > {primer3_outputfile_name}"
	result = subprocess.run(command,shell=True)

	#Collect output parameters from primer3 output file.
	with open(primer3_outputfile_name) as file:
		tm = []
		gc = []
		selfany = []
		selfend = []
		hairpin = []
		stability = []
		for line in file:
			key,value = line.strip().split('=')
			if (key=="PRIMER_LEFT_0_TM"):
				tm.append(float(value))
			elif (key=="PRIMER_LEFT_0_GC_PERCENT"):
				gc.append(float(value))
			elif(key=="PRIMER_LEFT_0_SELF_ANY_TH"):
				selfany.append(float(value))
			elif(key=="PRIMER_LEFT_0_SELF_END_TH"):
				selfend.append(float(value))
			elif(key=="PRIMER_LEFT_0_HAIRPIN_TH"):
				hairpin.append(float(value))
			elif(key=="PRIMER_LEFT_0_END_STABILITY"):
				stability.append(float(value))
	return [seqid,sourceseq,format(statistics.mean(tm),'.2f'),format(statistics.mean(gc),'.2f'),format(statistics.mean(selfany),'.2f'),format(statistics.mean(selfend),'.2f'),format(statistics.mean(hairpin),'.2f'),format(statistics.mean(stability),'.2f')]

#----------------------------------------------------------------------------------------------------

def run_primer3(consensus_list):
	return_list = []
	#Prepare primer3 input and output.
	primer3_inputfile_name = f"primer3_input_tmp.txt"
	primer3_outputfile_name = f"primer3_output_tmp.txt"

#	print(f"Primer_name\tPrimer_seq\tTm\tGC_content\tSelf_any\tSelf_end\tHairpin\t3\'Stability")

	for specname,start,end,consensus_sequence,consensus_sequence_up10,consensus_sequence_down10,start_pos,end_pos in consensus_list:
		seqid = specname+"_"+str(start)+"_"+str(end)
		seq = consensus_sequence
		upseq = consensus_sequence_up10
		downseq = consensus_sequence_down10

		bestcandidate_for = ["","",0,100,100,100,100,100]

		#Adding extra bases at most 5 nt to elongate the primer candidate sequence.
		for i in range(6):
			for j in range(6):
				upaddition = ""
				downaddition = ""
				if (i>0 and len(upseq)>=i):
					upaddition = upseq[len(upseq)-i:]
				if (j>0 and len(downseq)>=j):
					downaddition = downseq[:j]

				#Test forward sequence.
				candidate_seq = upaddition+seq+downaddition
				if (
					len(IUPACData.ambiguous_dna_values[candidate_seq[0]])>1 or
					len(IUPACData.ambiguous_dna_values[candidate_seq[-1]])>1
					):
					continue
				candidate_seqid = "f_"+seqid+"_"+str(i)+"_"+str(j)
				all_possible_sequence = expand_degenerate_bases(candidate_seq)
				results = run_primer3_seqlist(candidate_seqid,candidate_seq,all_possible_sequence,primer3_inputfile_name,primer3_outputfile_name)
				if (
					bestcandidate_for[0]=="" and
					float(results[2])>=55 and
					float(results[2])<=65 and
					float(results[4])<=50 and
					float(results[5])<=50 and
					float(results[6])<=50 and
					float(results[7])<=4
					):
					bestcandidate_for = results
				elif (
					float(results[2])>=55 and
					float(results[2])<=65 and
					float(results[4])<=float(bestcandidate_for[4]) and
					float(results[4])<=50 and
					float(results[5])<=float(bestcandidate_for[5]) and
					float(results[5])<=50 and
					float(results[6])<=float(bestcandidate_for[6]) and
					float(results[6])<=50 and
					float(results[7])<=float(bestcandidate_for[7]) and
					float(results[7])<=4
					):
					bestcandidate_for = results

		if (bestcandidate_for[0]!=""):
			return_list.append(bestcandidate_for)
#			print(*bestcandidate_for,sep="\t")
	cmd = 'rm primer3_input_tmp.txt primer3_output_tmp.txt'
	res = subprocess.run(cmd,shell=True)

	return return_list

#----------------------------------------------------------------------------------------------------

def import_mismatch(file_path):
    mismatch_list = [0,0,0,0,0]
    primer_name = ""
    with open(file_path, 'r') as file:
        counter = 0
        for line in file:
            columns = line.strip().split(',')
            if re.search(r'# Primer: .* 5\'',columns[0]):
                primer_name = re.findall(r'# Primer: (.*) 5\'',columns[0])[0]
            elif not(re.search(r'#',columns[0])):
                counter += 1
                mismatch_num = int(columns[4])+int(columns[5])
                if columns[6]==True:
                    mismatch_num += 1
                if mismatch_num > 3:
                    mismatch_num = 4
                mismatch_list[mismatch_num]+=1
    mismatch_list.append(primer_name)
    mismatch_list.append(counter)
    return mismatch_list

#----------------------------------------------------------------------------------------------------

def extractppresults(targetdirname,nontargetdirname):
	target_file_list = sorted(glob.glob(os.path.join(targetdirname,'*.txt')))
	nontarget_file_list = sorted(glob.glob(os.path.join(nontargetdirname,'*.txt')))
#	print("primer_name","Mismatch_0_fraction","Mismatch_0_1_fraction","Mismatch_0_2_fraction",sep='\t')
	return_list = []
	primer_list = []
	target_dict = {}
	nontarget_dict = {}
	for txtfile in target_file_list:
		mismatch_list = import_mismatch(txtfile)
		primer_list.append(mismatch_list[5])
		target_dict[mismatch_list[5]]=(
			round(mismatch_list[0]/mismatch_list[6],2),
			round((mismatch_list[0]+mismatch_list[1])/mismatch_list[6],2),
			round((mismatch_list[0]+mismatch_list[1]+mismatch_list[2])/mismatch_list[6],2)
			)
	for txtfile in nontarget_file_list:
		mismatch_list = import_mismatch(txtfile)
		nontarget_dict[mismatch_list[5]]=(
			round(mismatch_list[0]/mismatch_list[6],2),
			round((mismatch_list[0]+mismatch_list[1])/mismatch_list[6],2),
			round((mismatch_list[0]+mismatch_list[1]+mismatch_list[2])/mismatch_list[6],2)
		)
	for primer in primer_list:
		tmparray = (primer,target_dict[primer][0],target_dict[primer][1],target_dict[primer][2],nontarget_dict[primer][0],nontarget_dict[primer][1],nontarget_dict[primer][2])
		return_list.append(tmparray)

	return return_list

#----------------------------------------------------------------------------------------------------


if __name__ == '__main__':
	try:
		args = get_arguments()
		inputfile = args.input
		outputfolder = re.sub('$/','',args.outputfolder)
		infolder, infile = os.path.split(os.path.realpath(inputfile))
		outprefix = args.prefix
		taxid = args.taxid
		datadir = args.datadir
		taxkitoption = args.taxkitoption
		entthreshold = args.entropythreshold
		refseqid = args.baseposRefSeqid

		print('\n#00: Prepare directories.')
		cmd = 'mkdir -p {0}/01_sequence_splitting/tmp {0}/02_entropycalc {0}/03_primer_candidates'.format(outputfolder)
		print(cmd)
		res = subprocess.run(cmd,shell=True)

		############################
		# 01 Sequence Splitting    #
		############################
		print('\n#01: Sequence splitting.')
		print('\n#01-1: Prepare correspondence file of taxon ID and RefSeqID of the sequences in the input file.')
		#File check.
		nucl_gb_gz = datadir+"/nucl_gb.accession2taxid.gz"
		nucl_gb = datadir+"/nucl_gb.accession2taxid"
		nucl_gb_smaller = datadir+"/nucl_gb.accession2taxid_smaller_sorted.txt"
		if(not os.path.exists(nucl_gb_gz) and not os.path.exists(nucl_gb) and not os.path.exists(nucl_gb_smaller)):
			print('\nDownloading nucl_gb.accession2taxid.gz...')
			cmd = 'wget -P {0} https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz'.format(datadir)
			print(cmd)
			res = subprocess.run(cmd,shell=True)
		if(os.path.exists(nucl_gb_gz) and not os.path.exists(nucl_gb) and not os.path.exists(nucl_gb_smaller)):
			print('\nExtracting nucl_gb.accession2taxid.gz...')
			cmd = 'gunzip {0}/nucl_gb.accession2taxid.gz'.format(datadir)
			print(cmd)
			res = subprocess.run(cmd,shell=True)
		if(os.path.exists(nucl_gb) and not os.path.exists(nucl_gb_smaller)):
			print('\nMaking smaller and sorted version of nucl_gb.accession2taxid...')
			cmd = 'awk \'{{print $2, $3}}\' {0}/nucl_gb.accession2taxid | sort > {0}/nucl_gb.accession2taxid_smaller_sorted.txt'.format(datadir)
			print(cmd)
			res = subprocess.run(cmd,shell=True)
		if(os.path.exists(nucl_gb_smaller)):
			print('\nnucl_gb.accession2taxid_smaller_sorted.txt prepared.')

		print('\n#01-2: Make the list of Refseq ID of the sequences in the input file.')
		cmd = 'cat {0} | awk -F\' \' \'$0~/>/{{print $1}}\' | sed -e \'s/>//g\' | sort > {1}/01_sequence_splitting/tmp/input_refseqid_sorted.txt'.format(inputfile,outputfolder)
		print(cmd)
		res = subprocess.run(cmd,shell=True)
		cmd = 'join -1 1 -2 1 {0}/01_sequence_splitting/tmp/input_refseqid_sorted.txt {1}/nucl_gb.accession2taxid_smaller_sorted.txt -o 2.1,2.2 | sort -k2,2 > {0}/01_sequence_splitting/tmp/input_refseqid_taxid_sorted.txt'.format(outputfolder,datadir)
		print(cmd)
		res = subprocess.run(cmd,shell=True)

		print('\n#01-3: Prepare taxID of target organism group.')
		cmd = 'taxonkit list {0} --ids {1} --indent "" | sort > {2}/01_sequence_splitting/tmp/target_taxid_sorted.txt'.format(taxkitoption,taxid,outputfolder)
		print(cmd)
		res = subprocess.run(cmd,shell=True)
		cmd = 'join -1 1 -2 2 {0}/01_sequence_splitting/tmp/target_taxid_sorted.txt {0}/01_sequence_splitting/tmp/input_refseqid_taxid_sorted.txt -o 2.1 > {0}/01_sequence_splitting/tmp/target_refseqid.txt'.format(outputfolder)
		print(cmd)
		res = subprocess.run(cmd,shell=True)
		#Count number of target ids.
		idlist = outputfolder + "/01_sequence_splitting/tmp/target_refseqid.txt"
		line_count = 0
		with open(idlist,"r",encoding="utf-8") as f:
			line_count = sum(1 for _ in f)

		print('\n#01-4: Prepare target and non-target sequence fasta files.')
		cmd = 'seqkit grep -f {0}/01_sequence_splitting/tmp/target_refseqid.txt {1} -o {0}/01_sequence_splitting/{2}_target_mtDNA.fasta'.format(outputfolder,inputfile,outprefix)
		print(cmd)
		res = subprocess.run(cmd,shell=True)
		cmd = 'seqkit grep -v -f {0}/01_sequence_splitting/tmp/target_refseqid.txt {1} | seqkit sample -n {3} -o {0}/01_sequence_splitting/{2}_nontarget_mtDNA.fasta'.format(outputfolder,inputfile,outprefix,line_count)
		print(cmd)
		res = subprocess.run(cmd,shell=True)

		print('\n#01-5: Unify the start position and direction of mtDNA sequences.')
		target_fastafile_base = outprefix+"_target_mtDNA"
		if(os.path.exists(outputfolder+"/"+target_fastafile_base+"_setstart.fasta")):
			cmd = 'rm {0}/{1}_setstart.fasta'.format(outputfolder,target_fastafile_base)
			res = subprocess.run(cmd,shell=True)
		setstartposition(outputfolder,target_fastafile_base)
		cmd = 'mv {0}_target_mtDNA_setstart.fasta {1}/01_sequence_splitting/'.format(outprefix,outputfolder)
		res = subprocess.run(cmd,shell=True)

		print('\n#01-6: Align target sequence file using MAFFT.')
		cmd = 'mafft --thread -1 {0}/01_sequence_splitting/{1}_target_mtDNA_setstart.fasta > {0}/01_sequence_splitting/{1}_target_mtDNA_setstart_aligned.fasta'.format(outputfolder,outprefix)
		print(cmd)
		res = subprocess.run(cmd,shell=True)
		print("\nThis step would take a while if the target sequence file includes large number of entries...")

		############################
		# 02 Entropy calculation   #
		############################

		print('\n#02-1: Pick the first RefSeq ID of the representative of the target sequences.')
		posbase_refseqid = ""
		refseqidpath = outputfolder+"/01_sequence_splitting/tmp/target_refseqid.txt"
		if refseqid:
			posbase_refseqid = refseqid
		else:
			with open(refseqidpath,'r',encoding='utf-8') as file:
				refseqids = [line.strip() for line in file]
			if refseqids:
				posbase_refseqid = refseqids[0]
			else:
				print("\nError. Target sequence was not picked at all.")
				sys.exit()

		print('\n#02-2: Calculate entropy of each position on mtDNA.')
		alignedfile = outputfolder+"/01_sequence_splitting/"+outprefix+"_target_mtDNA_setstart_aligned.fasta"
		entropylistfile = outputfolder+"/02_entropycalc/"+outprefix+"_entropy_"+posbase_refseqid+".txt"
		entarray = entropy_calc(alignedfile, posbase_refseqid, entropylistfile, entthreshold)
		#File output
		enthigerthresholdfile = outputfolder+"/02_entropycalc/"+outprefix+"_entropy_"+posbase_refseqid+"_ent"+str(entthreshold)+".txt"
		with open(enthigerthresholdfile, 'w', encoding='utf-8') as f:
			f.write(f"Specname\tStart\tEnd\tStart_position_in_alignment\tEnd_position_in_alignment\tInformation_content\n")
			for array in entarray:
				line = '\t'.join(map(str,array))
				f.write(line+'\n')

		#############################
		# 03 Primer candidate pick  #
		#############################

		print('\n#03-1: Extract the best conserved 20 bp region from conserved regions.')
		primercandidatesfile = outputfolder+"/03_primer_candidates/"+outprefix+"_entropy_"+posbase_refseqid+"_ent"+str(entthreshold)+"_bestregions.txt"
		conservedregion_list = extractconservedregionseq(entarray,primercandidatesfile)

		print('\n#03-2: Extract consensus sequences from the conserved regions.')
		alignmentfile = outputfolder+"/01_sequence_splitting/"+outprefix+"_target_mtDNA_setstart_aligned.fasta"
		consensusseq_info = getconsensus_fromalignment(conservedregion_list, alignmentfile)
		#File output
		consensusseqfile = outputfolder+"/03_primer_candidates/"+outprefix+"_entropy_"+posbase_refseqid+"_ent"+str(entthreshold)+"_bestregions_seq.txt"
		with open(consensusseqfile, 'w', encoding='utf-8') as f:
			f.write(f"Specname\tStart\tEnd\tConsensus_sequence\tConsensus_sequence_up10\tConsensus_sequence_down10\tStart_position_in_alignment\tEnd_position_in_alignment\n")
			for array in consensusseq_info:
				line = str(array[3])+"\t"+str(array[0])+"\t"+str(array[1])+"\t"+str(array[2])+"\t"+str(array[4])+"\t"+str(array[5])
				f.write(line+'\n')

		print('\n#03-3: Run primer3 to design primer sequence candidates')
		primerinfo = run_primer3(consensusseq_info)
		#File output
		primeroutfile = outputfolder+"/03_primer_candidates/"+outprefix+"_entropy_"+posbase_refseqid+"_ent"+str(entthreshold)+"_bestregions_seq_primer3.txt"
		with open(primeroutfile, 'w', encoding='utf-8') as f:
			f.write(f"Primer_name\tPrimer_seq\tTm\tGC_content\tSelf_any\tSelf_end\tHairpin\t3\'Stability\n")
			for array in primerinfo:
				line = '\t'.join(map(str,array))
				f.write(line+'\n')

		print('\n#03-4: Run PrimerProspector for selectivity filtering.')
		ppinputfile = "03_primer_candidates/"+outprefix+"_entropy_"+posbase_refseqid+"_ent"+str(entthreshold)+"_bestregions_seq_primer3_ppinput.txt"
		count = 0
		with open(ppinputfile, 'w', encoding='utf-8') as f:
			for array in primerinfo:
				count = count + 1
				zeropaddedcount = f"{count:03}"
				primercandseq = str(array[1])
				line = zeropaddedcount+array[0]+"\t"+primercandseq
				f.write(line+'\n')
		targetseqfile = "01_sequence_splitting/"+outprefix+"_target_mtDNA.fasta"
		nontargetseqfile = "01_sequence_splitting/"+outprefix+"_nontarget_mtDNA.fasta"

		#Target specificity
		cmd = '. $CONDA_PREFIX/../../etc/profile.d/conda.sh && conda activate {0} && analyze_primers.py -f {1} -P {2} -o 03_primer_candidates/{3}_primerprospector_target/ && conda deactivate'.format(args.ppenv,targetseqfile,ppinputfile,outprefix)
		print(cmd)
		res = subprocess.run(cmd,shell=True)
		#Non-target specficity
		cmd = '. $CONDA_PREFIX/../../etc/profile.d/conda.sh && conda activate {0} && analyze_primers.py -f {1} -P {2} -o 03_primer_candidates/{3}_primerprospector_nontarget/ && conda deactivate'.format(args.ppenv,nontargetseqfile,ppinputfile,outprefix)
		print(cmd)
		res = subprocess.run(cmd,shell=True)

		mismatch_list = extractppresults(outputfolder+"03_primer_candidates/"+outprefix+"_primerprospector_target",
										outputfolder+"03_primer_candidates/"+outprefix+"_primerprospector_nontarget")

		ppoutfile = outputfolder+"/03_primer_candidates/"+outprefix+"_entropy_"+posbase_refseqid+"_ent"+str(entthreshold)+"_bestregions_seq_primer3_ppout.txt"
		with open(ppoutfile, 'w', encoding='utf-8') as f:
			for array in mismatch_list:
				line = '\t'.join(map(str,array))
				f.write(line+'\n')

		ppoutselectedfile = outputfolder+"/03_primer_candidates/"+outprefix+"_entropy_"+posbase_refseqid+"_ent"+str(entthreshold)+"_bestregions_seq_primer3_ppout_candidates.txt"
		with open(ppoutfile, 'w', encoding='utf-8') as f:
			for array in mismatch_list:
				if (array[2]>=0.85 and array[5]<=0.15):
					line = '\t'.join(map(str,array))
					f.write(line+'\n')

	except KeyboardInterrupt:
		sys.exit()