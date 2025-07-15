import sys
import re
import math
from Bio import AlignIO

#Output as follows:
#Entry	startpos_in_entry	endpos_in_entry	corresponding_startpos_in_alignment	corresponding_endpos_in_alignment	average_entropy


args = sys.argv
if len(args)!=3:
	print('Usage: args[0] multiple_alignment_file outputentry_id')
else:
	file = args[1]
	outputname = args[2]
	align = AlignIO.read(file,"fasta")
	alignmentlength = align.get_alignment_length()
	samplenum = 0
	output_ind = 0
	for i in range(len(align)):
		if (re.match(r'Consensus',align[i].id)==None):
			samplenum = samplenum + 1
		if (re.match(outputname,align[i].id)):
			output_ind = i
#	print(output_ind)
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

	tmppos = 0
	while tmppos < poscount-20:
		tt = 0
		for j in range(tmppos,tmppos+20):
			tt = tt + entlist[j]
		tt = tt/20
		print(align[output_ind].id,tmppos+1,tmppos+20,alignmentpos[tmppos+1],alignmentpos[tmppos+20],tt,sep='\t')
		tmppos = tmppos + 1