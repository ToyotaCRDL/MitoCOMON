from Bio import SeqIO
import os
import sys
import re
import subprocess
from subprocess import PIPE

args = sys.argv
if len(args)!=3:
	print('Usage: setstartposition_mtDNA.py contig_before_rotation.fasta prefix')
	sys.exit()

fastafile = args[1]
prefix = args[2]

startpos = 0
endpos = 0
strand = 0
score = 0

n = 0

fasta = SeqIO.parse(fastafile,'fasta')
for r in fasta:
	n = n+1
	print("Started analyzing contig",n)
	seq = r.seq
	seqid = r.id
#	pattern = re.compile(r'.*.circularised')
#	match = pattern.search(seqid)
#	if not match:
#		print(f"Contig {r.id} was not circularised. Skip.")
#		continue

	with open('./tmptrnascantarget.fasta','w') as f:
		f.write(r.format('fasta-2line'))

	print("Running tRNAscan-SE to detect tRNAs...")
	cmd = 'tRNAscan-SE {} -M vert'.format(os.path.abspath('./tmptrnascantarget.fasta'))
	res = subprocess.Popen(cmd,shell=True,stdout=PIPE)

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

	outfilename = prefix+".fasta"

	if strand == 0:
		#Rotate and reanalyze the sequence for the case that the contig start in the middle of tRNA-Phe.
		print("Retry by rotating the contig.")
		seq = seq[1000:len(r)]+seq[0:1000]
		reanalyzeseq = ">rotated_entry\n"+str(seq)
		with open('./tmptrnascantarget.fasta','w') as f:
			f.write(reanalyzeseq)
		print("Running tRNAscan-SE to detect tRNAs...")
		cmd = 'tRNAscan-SE {} -M vert'.format(os.path.abspath('./tmptrnascantarget.fasta'))
		res = subprocess.Popen(cmd,shell=True,stdout=PIPE)
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
		if strand == 0:
			print("Couldn't find tRNA-Phe in contig",n)
			continue

	print("tRNA-Phe position information:")
	print("contig",n,startpos,endpos,strand)

	rotated_seq = ''
	if strand == 1:
		subseq1 = seq[0:startpos-1]
		subseq2 = seq[startpos-1:len(r)]
		rotated_seq = subseq2+subseq1
	elif strand == -1:
		subseq1 = seq[0:startpos]
		subseq2 = seq[startpos:len(r)]
		rotated_seq = subseq1.reverse_complement()+subseq2.reverse_complement()


	idname = ">contig"+str(n)
	
	if n==1:
		with open(outfilename,'w') as f:
			f.write(idname)
			f.write("\n")
			f.write(str(rotated_seq))
			f.write("\n")
	else:
		with open(outfilename,'a') as f:
			f.write(idname)
			f.write("\n")
			f.write(str(rotated_seq))
			f.write("\n")


	print("Wrote the rotated contig sequence in ",outfilename,".",sep='')


