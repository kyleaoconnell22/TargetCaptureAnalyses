import sys
from Bio import AlignIO
from Bio.Align import AlignInfo
import os
import subprocess as sp
import shutil
import argparse
import random


'''
going to take alignment and create consensus from phased alleles
Input is a phylip aln and also a list of samples
make sure the aln has 
sample1_0 ATGTACA
sample1_1 ATGTACA
sample2_0 ATGTACA
sample2_1 ATGTACA
 ... 

The sample file should just be: 
sample1
sample2
 ...
The script will only keep the samples in this file for the consensus, so you can also use this file to filter out unwanted individuals.

#The final output here is a fasta file in gblocks format. I used Geneious to conver to Phylip for RAxML but you could add a module with Bio Python to convert to other formats.
Kyle O'Connell
kyleaoconnell22@gmail.com
March 2020
'''


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--in_dir", required=True, help="REQUIRED: The full path to a directory which contains the input alignment files.")
    parser.add_argument("-s", "--in_sp", required=True, help="REQUIRED: list of all samples you want to include in the consensus aln")
    parser.add_argument("-p", "--in_phy", required=True, help="REQUIRED: input phylip aln")
    return parser.parse_args()

def in_samples(in_dir,in_sp):
	
	#open list of samples without the allele designation
	fh_sp = open(in_sp,'r')
	sp_list = []
	for line in fh_sp:
		line=line.strip()
		sp_list.append(line)
	
	sp_set = set(sp_list)
	return sp_set
	
def write_out(in_phy):
	out_file = in_phy.split('.')[0]+'_consensus.fasta'
	return out_file
	
def split_alleles(in_dir,in_phy,in_sp):
	os.chdir(in_dir)

	#open phylip aln should already be aligned and phased alleles
	fh_aln = open(in_phy, 'r')

	#call the list of samples from above
	samples = in_samples(in_dir,in_sp)
	
	seqs_0 = [] #create empty dicts to write all 0 alleles 
	seqs_1 = [] #empty dict for all 1 alleles
	#write seqs to dicts
	for seq in fh_aln:
		if seq.split(' ')[0].endswith('0'):
			seq = seq.strip()
			seqs_0.append(seq)
		elif seq.split(' ')[0].endswith('1'):
			seq = seq.strip()
			seqs_1.append(seq)
	
	#sort our seq lists, I did this in dicts but had issues with iterating below, 
	#so I switched back to lists
	#could be an issue with memory if the alignment is super huge, but mine has 400k bp
	
	s0 = []
	s1 = []
	
	seqs_0.sort()
	seqs_1.sort()
	###check for missing alleles
	for i,j in zip(seqs_0,seqs_1):
		s0.append(i.split(' ')[0].split('_')[0])
		s1.append(j.split(' ')[0].split('_')[0])
		print i.split(' ')[0], j.split(' ')[0]
	for i,j in zip(seqs_0,seqs_1):
		if i.split(' ')[0].split('_')[0] not in s1:
			print i.split(' ')[0].split('_')[0]
		elif j.split(' ')[0].split('_')[0] not in s0:
			print j.split(' ')[0].split('_')[0]
	
	#sort seq lists by sample name, can add module to check that they are the same, but I have already done that
	#need to be in the same order with the same number of 0 and 1 samples
	#iterate through the two lists
	for i,j in zip(seqs_0,seqs_1):
		#write temp fasta file with each seq
		temp_out = i.split(' ')[0].split('_')[0]+'.fasta'
		fh_temp = open(temp_out,'a')
		fh_temp.write('>'+i.split(' ')[0]+'\n'+i.split(' ')[1]+'\n''>'+j.split(' ')[0]+'\n'+j.split(' ')[1])

def call_consensus(in_dir,in_phy):
	os.chdir(in_dir)
	
	###open outfile####
	outfile = in_phy.split('.')[0]+'_consensus.fasta.txt'
	fh_out = open(outfile,'a')
		
	for filetype in os.listdir('.'):
		if filetype.endswith('.fasta'):
			alignment = AlignIO.read(filetype, 'fasta')
			summary_align = AlignInfo.SummaryInfo(alignment)
			print 'generating consensus of', filetype
			consensus = summary_align.gap_consensus()
			fh_out.write('>'+filetype.split('.')[0]+'\n'+str(consensus)+'\n')
			os.remove(filetype)
		elif filetype.endswith('.fasta.txt'):
			newname = filetype.split('.')[0]+'.fasta'
			os.rename(filetype,newname)


def remove_fails(in_dir,in_phy):
	os.chdir(in_dir)
	
	new_out = in_phy.split('.')[0]+'_consensus_trimmed.fasta'
	fh_out=open(new_out,'a')
	###open outfile####
	outfile = write_out(in_phy)
	alignment = AlignIO.read(outfile, 'fasta')
	
	
	for item in alignment:
		seq = str(item.seq)
		seq=seq.replace('X','-')
		dash = item.seq.count('-')
		leng = len(item.seq)
		perc = float(dash)/(len(item.seq))
		#round to two decimal places
		perc = round(perc,2)
		if perc < 0.75:
			fh_out.write(str('>'+item.id+'\n'+seq+'\n'))
		else:
			print item.id
			
				
def trim_gblocks(in_dir,in_phy):
	os.chdir(in_dir)
	
	in_aln = in_phy.split('.')[0]+'_consensus_trimmed.fasta'
	
	call_string = "gblocks {0} -b2=.1 -b3=8 -b4=10 -b5=a -t=d".format(in_aln)
	proc = sp.call(call_string, shell=True)
	
####at this point I just converted the output fasta to phy using Genious. 
#but you could deal with the interleaved output and convert to phy if I want to add another module.


def main():
	#define the arguments
	args = get_args()
	in_samples(args.in_dir,args.in_sp)
	write_out(args.in_phy)
	split_alleles(args.in_dir,args.in_phy,args.in_sp)
	call_consensus(args.in_dir,args.in_phy)
	remove_fails(args.in_dir,args.in_phy)
	trim_gblocks(args.in_dir,args.in_phy)
	
if __name__ == '__main__':
    main()	
