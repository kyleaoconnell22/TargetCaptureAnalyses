#!/usr/bin/python

import argparse
import os
import subprocess as sp
import shutil
from datetime import datetime
from Bio import SeqIO
from Bio.Seq import Seq
import numpy
import numpy as np


'''
Lots of this script was adapted from SuperCrunch by portik@github

This assumes the data have been phased, but there is also an unphased version in this repo

Notice that
j = number of '-' char
i = total len of aln

Kyle O'Connell
kyleaoconnell22@gmail.com
3/20/19

python capture_trim_summary_v2.py -i . -f fasta -a all
'''
#define the user inputs 
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--in_dir", required=True, help="REQUIRED: The full path to a directory which contains the input alignment files.")
    parser.add_argument("-f", "--out_format", required=True, choices=["fasta","nexus","phylip"], help="REQUIRED: Specify the output file format for trimmed alignment.")
    parser.add_argument("-a", "--analysis", required=True, choices=["trimal", "gblocks","all_trimming","summary", "all"], help="REQUIRED: Specify the desired trim method and or summary function.")
    return parser.parse_args()

#calc averages
def get_average(array):
	ave = np.average(array)
	return round(ave,2)

#calc percentage for missing data calc
def percentage(part, whole):
	if whole == 0:
		pass
	else:
  		return round(100 * float(part)/float(whole),2)

#remove indiv with only a few bp of seq or all missing data
#this is especially a problem with the gblocks trimming where the whole seq gets trimmed and it leaves just dashes
def remove_fails(file,suffix):
	with open(file, "rU") as handle:
		prefix = file.split('.')[0]
		prefix = prefix.split('_')[0] + '_' + prefix.split('_')[1] + '_' + prefix.split('_')[2]
		#take the temp files from the trimming below and write the filtered data to a new file with the 
		#suffix either trimal or gblocks, depending on the func call (suffix)
		outname = prefix + '_' + suffix + '.fasta'
		fh_out = open(outname, 'a')
		#parse the fasta files output by the trimming funcs below
		for record in SeqIO.parse(handle, "fasta"):
			#start counts to calc the proportion of '-'
			j = 0
			i = 0
			
			#if dash char count to j
			for char in record.seq:
				if char == '-':
					j = j + 1
			i = i + len(record.seq) #calc len of sequence
			missing = percentage(j,i) #calc perc missing (prop dash char)
			if missing > 80.0: #if we have more than 80% missing, throw this sample out
				pass
			else: #if at least 21% complete (in my data all the samples except the ones that really failed where at least 50%)
				seq = str(record.seq) #write info to output file
				fh_out.write('>' + record.id + '\n')
				fh_out.write(seq + '\n')
		os.remove(file)
			

#replace 'n' with '-'
'''
this and the next function replace the masked based in the alns from secapr, even with another pipeline
the script expects you to run this so its best to just leave it, even though it takes up more space in the end
'''
def replace_dash(file,out):
	#open sequence file
	with open(file, "rU") as handle:
		#parse using SeqIO
		for record in SeqIO.parse(handle, "fasta"):
			#grab the name without the fasta
			#write the fasta line to the outfile
			out.write('>' + record.id + '\n')
			seq = str(record.seq) #make seq.id a str
			seq = seq.upper() #make the bp upper case
			out.write((seq.replace('N','-')) + '\n') #replace the masked bases (Ns) and replace with - so they get trimmed below
	
#replaced masked based from secapr with '-' so that trimal will filter failed samples out
def replace_Ns(dir):
	#cd into unreplaced fasta dir
	os.chdir(dir)
	
	#check to see if ns have already been replaced because otherwise it crashes the script
	if os.path.exists("Output_replaced"):
		pass
	else:
	#if ns have not been replaced, keep going
	#iterate through fastas 
		for filetype in os.listdir('.'):
			#only open fasta files
			if filetype.endswith(".fasta"):
				#grab name for output file
				outfile = filetype.split('.')[0] + '_' + 'replaced' + '.fasta'
				#open output file for appending
				fh_out = open(outfile, 'a')
				#call the function above replace_dash
				replace_dash(filetype,fh_out)

#move all replaced files to a new directory, where all the rest of the script will happen from
def move_replaced(dir):
	os.chdir(dir)
	out_dir = "Output_replaced"
	if not os.path.exists(out_dir):
		os.mkdir(out_dir)
    #iterate through files in current directory   
	for filetype in os.listdir('.'):
		#only move replaced files to the output dir
		if filetype.endswith("_replaced.fasta"):
			shutil.move(filetype, out_dir)
	
#trimal module
'''
you can edit this with different params, such replacing -gt 0.1 to -gt 0.5 or -gappyout
but I found in my alns it makes very little difference
'''		
def trim_align_trimal(dir, out_format):
	print 'Starting trimal routine:'
	#make sure you are in the Output_replaced dir
	os.chdir(dir)
	for filetype in os.listdir('.'):
		print "\n\nTrimming alignment {} with trimal".format(filetype)
		prefix = filetype.split('.')[0]
		#this is where it can handle different types of inputs, although I am only using fasta right now
		if out_format == "fasta":
			call_string = "trimal -in {0} -out {1}_temp.fasta -fasta -gt 0.05".format(filetype, prefix)
		elif out_format == "nexus":
			call_string = "trimal -in {0} -out {1}_temp.nex -nexus -gt 0.1".format(filetype, prefix)
		elif out_format == "phylip":
			call_string = "trimal -in {0} -out {1}_temp.phy -phylip_paml -gt 0.1".format(filetype, prefix)
		#call it!
		proc = sp.call(call_string, shell=True)
	
	#call the remove fails function on these trimmed seqs
	'''
	I found that several of the seqs would trim, and would leave like 50 bp or have almost nothing, so here it calls
	the func above to remove these seqs
	'''
	for filetype in os.listdir('.'):
		if filetype.endswith('_temp.fasta'):
			remove_fails(filetype, 'trimal')
		
	#move output section
	#this could probably just be one function, but I would need to edit more things for that
	out_dir = "Output_trimal"
	if not os.path.exists(out_dir):
		os.mkdir(out_dir)
	
	for filetype in os.listdir('.'):
		#only move remove failed (filtered) files to the output dir
		if filetype.endswith("_trimal.fasta") or filetype.endswith("_trimal.nexus") or filetype.endswith("_trimal.phylip_paml"):
			shutil.move(filetype, out_dir)	
	for output in os.listdir('.'):
		if output == 'Output_trimal':
			shutil.move(output, '..')
	#have to cd out in order for gblocks to work
	#this is only an issue with 'all' or 'all_trimming'
	os.chdir('..')
	

#Gblocks Trimming
def trim_gblocks(dir,format):
	print 'Starting gblocks routine:'
	os.chdir(dir)
	for file in os.listdir('.'):
		if format == 'fasta':
			print "\n\nTrimming alignment {} with gblocks".format(file)
			call_string = "gblocks {0} -b2=.1 -b3=8 -b4=10 -b5=a -t=d".format(file)
			proc = sp.call(call_string, shell=True)
		else:
			print "gblocks can't trim other input formats"
	
	for filetype in os.listdir('.'):
		if filetype.endswith('.fasta-gb'):
			remove_fails(filetype, 'gblocks')
	
	#move outputs
	out_dir = "Output_gblocks"
	if not os.path.exists(out_dir):
		os.mkdir(out_dir)
	
	for filetype in os.listdir('.'):
		if filetype.endswith('.fasta-gb.htm') or filetype.endswith('fasta-gb'):
			os.remove(filetype)
		
	for filetype in os.listdir('.'):
		#only move replaced files to the output dir
		if filetype.endswith("_gblocks.fasta"):
			shutil.move(filetype, out_dir)	
	for x in os.listdir('.'):
		if x == 'Output_gblocks':
			shutil.move(x, '..')
	os.chdir('..')

#Summary routine
def summarize_alignments(out_dir,suffix):
	#write output with trimal or gblocks summary
	#it will run this function on both output dirs
	sum_out = suffix + '_summary.txt'
	fh_out = open(sum_out, 'a')
	
	if os.path.exists(out_dir):
		os.chdir(out_dir)
		# task 1 : count files in dir to know how many loci are left after trimming
		filenum = [] #list to count number of files
		indnum = [] #count num of ind per aln
		i = 0 #count total bp of sequence captured
		j = 0 #count '-' to find % missing data across alns
		seqs = [] #count length of seqs
		namelist = [] #count all names in all files
		namecounts = [] #count unique names in files
		
		for filetype in os.listdir('.'):
			#count number of alignments = #loci captured across samples
			filenum.append(filetype)
					
			#use SeqIO to read in the fasta files easier
			with open(filetype, "rU") as handle:
				seq_temp = [] #used to calc average length of each aln, rather than of each seq
				for record in SeqIO.parse(handle, "fasta"):
					if record.id.startswith("_"):
						pass
					else:
						namelist.append(record.id) #write all names to namelist
										
						i = i + len(record.seq) #calc total bp in dataset
						seq_temp.append(len(record.seq)) #write just the len of the aln, not every seq
					
					for char in record.seq:
						if char == '-': 
							j = j + 1 #calc the number of dashes, we will use this to find missing data below
					
				ave = np.average(seq_temp) #this is average of all seqs in this aln
				#write to seqs list above just one len for the whole locus so each locus has one entry in the final list
				seqs.append(ave) 
				
		
		#iterate through the namelist		
		for name in sorted(namelist):
			#write the full counts to another list
			namestats = name + '\t' + str(namelist.count(name))
			namecounts.append(namestats)
		
		##### Write the summary output stats #####		
		#write number of loci retained to output
		print 'there are {0} files in the {1} directory'.format(len(filenum),suffix)
		fh_out.write('there are {0} files in the {1} directory'.format(len(filenum),suffix))
		fh_out.write('\n\n')
		
		#write total bp in final alns
		fh_out.write('the alignments have {0} total bps'.format(i))
		fh_out.write('\n\n')
		
		#write % missing data ('-' char) in final alns
		fh_out.write('the alignments have {0}% missing data across all seqs'.format(percentage(j,i)))
		fh_out.write('\n\n')
		
		#write ave len of each alignment
		fh_out.write('The average len of a locus = {0} bp'.format(get_average(seqs[1:])))
		fh_out.write('\n\n')
		fh_out.write('The range of lengths = {0} to {1} bp'.format(min(seqs),max(seqs)))
		fh_out.write('\n\n')
		
		#get unique values of counts (using set and then converting back to list for sorting)
		nameset = set(namecounts)
		#convert back to list to sort it
		namecountsuniq = sorted(list(nameset))
		fh_out.write('{0} samples are represented in these alignments'.format(len(namecountsuniq)))
		fh_out.write('\n\n')
		
		#write average individual per alignment
		fh_out.write('average number of samples in each alignment = {0}'.format(len(namelist)/len(filenum)))
		fh_out.write('\n\n')
		
		#write header for indiv coverage list
		fh_out.write('sample_name' + '\t' + 'captured_loci' + '\n')
		
		#iterate through and write the capture success per sample to the output file
		#indiv 1 \t #loci present
		for item in namecountsuniq:
			fh_out.write(item + '\n')
			
	#if the dir doesn't exist then don't do anything
	#for example if you only trimmed with trimal not gblocks, then it will skip the gblocks summary routine
	else:
		pass
	os.chdir('..')

def main():
	#define the arguments
	args = get_args()
	#change into the base dir with the secapr output alignments
	os.chdir(args.in_dir)
	#replace Ns and move to output direct comment this line out if not needed
	#if commented out make sure to replace the line below that says Output_replaced with 'args.in_dir'
	replace_Ns(args.in_dir)
	move_replaced(args.in_dir)
	
	#different analysis options
	if args.analysis == "gblocks":
		trim_gblocks("Output_replaced",args.out_format)
	elif args.analysis == "trimal":
		trim_align_trimal("Output_replaced",args.out_format)
	elif args.analysis == "all_trimming":
		trim_gblocks("Output_replaced",args.out_format)
		trim_align_trimal("Output_replaced",args.out_format)
	elif args.analysis == "summary":
		summarize_alignments("Output_trimal",'trimal')
		summarize_alignments("Output_gblocks",'gblocks')
	elif args.analysis == "all":
		trim_gblocks("Output_replaced",args.out_format)
		trim_align_trimal("Output_replaced",args.out_format)
		summarize_alignments("Output_trimal",'trimal')
		summarize_alignments("Output_gblocks",'gblocks')

if __name__ == '__main__':
    main()
