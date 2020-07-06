#!/usr/bin/env python2
def depth_masker(infile,outfile,min_depth,max_depth,map_qual):
	yui=open(infile,'r')
	new=open(outfile,'a')
	chromo=''
	twritten=0
	written=0
	total=0
	ok=0
	passed=0
	tpassed=0
	skip=0
	temporary=''
	last_pos=0
	for line in yui:
		if line[0]=='#':
			continue
		total+=1
		if total%50000000==0:
			print total, 'sites analyzed'
		split=line.split()
		if chromo!=split[0]:
			if chromo!='':
				start=0
				for full in range(len(temporary)/60):
					new.write(temporary[start:start+60]+'\n')
					start+=60
				if len(temporary)%60!=0:
					new.write(temporary[start:]+'\n')
				written+=len(temporary)
				temporary=''
			twritten+=written
			tpassed+=passed
			passed=0
			written=0
			skip=0
			last_pos=0
			new.write('>'+split[0]+'\n')
		while last_pos+1<int(split[1]):
			last_pos+=1
			temporary+='N'
		if skip>0:
			skip-=1
			chromo=split[0]
			last_pos=int(split[1])
			continue
		split2=split[7].split(';')
		pos_dict=dict()
		ok=0
		if len(temporary)>65:
			new.write(temporary[:60]+'\n')
			written+=60
			temporary=temporary[60:]
		for info in split2:
			info_split=info.split('=')
			try:
				pos_dict[info_split[0]]=info_split[1]
			except IndexError:
				continue
		try:
			(depth,MQ,AF,AC)=(pos_dict['DP'],pos_dict['MQ'],pos_dict['AF1'],pos_dict['AC1'])
		except KeyError:
			temporary+='N'
			chromo=split[0]
			continue
		if (min_depth<=int(depth)<=max_depth) and (float(MQ)>=map_qual):
			if 'INDEL' in split[7]:
				if (',' not in split[4]) and (float(AF)>0.5):
					temporary=temporary[:-1]+split[4]
				else:
					temporary=temporary[:-1]+split[3]
				skip=len(split[3])-1
				continue
			elif AC=='0':
				passed+=1
				temporary+=split[3]
				ok=1
			elif (AC=='1') and (len(split[4])==1):
				if (float(AF)>0.5) and (split[4]!='.'):
					temporary+=split[4]
				else:
					temporary+=split[3]
				passed+=1
				ok=1
		if ok==0:
			temporary+='N'
		chromo=split[0]
		last_pos=int(split[1])
	start=0
	for full in range(len(temporary)/60):
		new.write(temporary[start:start+60]+'\n')
		start+=60
	if len(temporary)%60!=0:
		new.write(temporary[start:]+'\n')
	written+=len(temporary)
	twritten+=written
	tpassed+=passed
	print total,'snp positions in the file'
	print twritten,'positions written'
	print tpassed,'positions satisfied the criterion'
	yui.close()

	
from optparse import OptionParser
import os,sys
parser = OptionParser()
parser.add_option('-v','--vcf_file', dest='vcf_file',help='Absolute path to vcf file to be converted [required]',type='str',default='')
parser.add_option('-o','--outfile', dest='outfile',help='Absolute path to the output fasta file [required]',type='str',default='')
parser.add_option('-d','--min_depth', dest='min_depth',help='Sites with lower depth are written as "N"',type='int',default=0)
parser.add_option('-D','--max_depth', dest='max_depth',help='Sites with higher depth are written as "N"',type='int',default=0)
parser.add_option('-Q','--map_quality', dest='map_quality',help='Sites with lower mapping quality are written as "N"',type='int',default=0)

(options,args)=parser.parse_args()
if (options.vcf_file=='') or (options.outfile==''):
	print '\ndepth_masker\nRequired filed(s) not supplied\n# Converts vcf file to consensus fasta.\n# The script can be used freely and comes with no guarantee whatsoever\n'
	parser.print_help()
	sys.exit(1)

import time
start_time=time.time()
depth_masker(options.vcf_file,options.outfile,options.min_depth,options.max_depth,options.map_quality)
end_time=time.time()
total_time=int(end_time-start_time)
print 'Analyses finished in',str(total_time/(60*60))+'hr',str((total_time%(60*60))/60)+'min',str((total_time%(60*60))%60)+'sec'
