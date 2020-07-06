#!/usr/bin/env python
def paired_quality_filter(infile1, infile2, paired1,paired2,unpaired1,unpaired2,quality,length):
	yui1=open(infile1,'r')
	yui2=open(infile2,'r')
	newp1=open(paired1,'a')
	newp2=open(paired2,'a')
	newu1=open(unpaired1,'a')
	newu2=open(unpaired2,'a')
	end=0
	total=0
	total_paired=0
	unpaired1=0
	unpaired2=0
	while end==0:
		head1=yui1.readline()
		if head1=='':
			end=1
			break
		total+=1
		if total%200000==0:
			print total,'reads analyzed...'
		seq1=yui1.readline()
		strand1=yui1.readline()
		qual1=yui1.readline()
		qual_seq=''
		qpos=0
		qual_list1=[]
		start=0
		on=0
		while qpos<len(qual1)-1:
			if ord(qual1[qpos])-33<quality:
				if on==1:
					qual_list1.append([qpos-start,start])
				on=0
			elif on==0:
				on=1
				start=qpos
				
			qpos+=1
		if on==1:
			qual_list1.append([qpos-start,start])
		qual_list1.sort(reverse=True)
		head2=yui2.readline()
		if head1[:-20]!=head2[:-20]:
			print 'Reads not properly paired!'
			end=1
			break
		seq2=yui2.readline()
		strand2=yui2.readline()
		qual2=yui2.readline()
		qpos=0
		qual_list2=[]
		start=0
		on=0
		while qpos<len(qual2)-1:
			if ord(qual2[qpos])-33<quality:
				if on==1:
					qual_list2.append([qpos-start,start])
				on=0
			elif on==0:
				on=1
				start=qpos
				
			qpos+=1
		if on==1:
			qual_list2.append([qpos-start,start])
		qual_list2.sort(reverse=True)
		if ((len(qual_list1)>0) and (qual_list1[0][0]>=length)) and ((len(qual_list2)>0) and (qual_list2[0][0]>=length)):
			total_paired+=1
			newp1.write(head1+seq1[qual_list1[0][1]:qual_list1[0][1]+qual_list1[0][0]]+'\n'+strand1+qual1[qual_list1[0][1]:qual_list1[0][1]+qual_list1[0][0]]+'\n')
			newp2.write(head2+seq2[qual_list2[0][1]:qual_list2[0][1]+qual_list2[0][0]]+'\n'+strand2+qual2[qual_list2[0][1]:qual_list2[0][1]+qual_list2[0][0]]+'\n')
		elif (len(qual_list1)>0) and (qual_list1[0][0]>=length):
			unpaired1+=1
			newu1.write(head1+seq1[qual_list1[0][1]:qual_list1[0][1]+qual_list1[0][0]]+'\n'+strand1+qual1[qual_list1[0][1]:qual_list1[0][1]+qual_list1[0][0]]+'\n')
		elif (len(qual_list2)>0) and (qual_list2[0][0]>=length):
			unpaired2+=1
			newu2.write(head2+seq2[qual_list2[0][1]:qual_list2[0][1]+qual_list2[0][0]]+'\n'+strand2+qual2[qual_list2[0][1]:qual_list2[0][1]+qual_list2[0][0]]+'\n')
	yui1.close()
	yui2.close()
	newp1.close()
	newp2.close()
	newu1.close()
	newu2.close()
	print 'Total input paired reads:',total
	print 'Total filtered paired reads:',total_paired
	print 'Filtered unpaired read1:',unpaired1
	print 'Filtered unpaired read2:',unpaired2

from optparse import OptionParser
import os,sys
parser = OptionParser()
parser.add_option('-i','--infile1', dest='infile1',help='Absolute path to the input fastq1 file [required]',type='str',default='')
parser.add_option('-j','--infile2', dest='infile2',help='Absolute path to the input fastq2 file [required]',type='str',default='')
parser.add_option('-1','--paired1', dest='paired1',help='Absolute path to the filtered paired1 fastq file [required]',type='str',default='')
parser.add_option('-2','--paired2', dest='paired2',help='Absolute path to the filtered paired2 fastq file [required]',type='str',default='')
parser.add_option('-u','--unpaired1', dest='unpaired1',help='Absolute path to the filtered unpaired1 fastq file [required]',type='str',default='')
parser.add_option('-v','--unpaired2', dest='unpaired2',help='Absolute path to the filtered unpaired2 fastq file [required]',type='str',default='')
parser.add_option('-q','--quality', dest='quality',help='Minimum phred quality to extract [Default: 20]',type='int',default=20)
parser.add_option('-l','--length', dest='length',help='Minimum length to keep [Default: 50]',type='int',default=50)

(options,args)=parser.parse_args()
if (options.infile1=='') or (options.infile2=='') or (options.paired1=='') or (options.paired2=='') or (options.unpaired1=='') or (options.unpaired2==''):
	print '\nRequired filed(s) not supplied\n The script can be used freely and comes with no guarantee whatsoever.\n##The program filters fastq file using quality and length information, and then extract the paired reads.\n## Unpaired reads that passed the filtering thresholds are also written\n'
	parser.print_help()
	sys.exit(1)
paired_quality_filter(options.infile1,options.infile2,options.paired1,options.paired2,options.unpaired1,options.unpaired2,options.quality,options.length)
