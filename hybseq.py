#!/usr/bin/env python

import os
import sys
import subprocess
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt



import argparse


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("-s", "--skim", type=str, help="skimming data", required=True)
	parser.add_argument("-t", "--trans", type=str, help="transcriptome data", required=True)
	parser.add_argument("-e", "--exon_prop", type=float, help="proportion of exon", default=0.7)
	args = parser.parse_args()
	skim_contig = args.skim
	transciptome = args.trans
	exon = args.exon_prop



print "***Starting of the hybseq pipeline***"


def cdhit(infile_cdhit, outfile_cdhit, cdhit_out_dir, c, tag):
	"""Clustering sequences with 80% or greater similarity"""
	c= str(c)
	cdhit_cmd = "cd-hit-est -i " + infile_cdhit + " -o " + outfile_cdhit + " -d 0 -c 1.0 -n 8 -p 1"
	screen_cdhit = tag + "_cdhit_" + c + "_screen"

	try:
		with open((cdhit_out_dir + screen_cdhit), 'w') as cdhit_screen:
			p = subprocess.Popen(cdhit_cmd,	stdout=cdhit_screen, shell=True)
			p.communicate()
		print ("\t\tcdhit completed.")
	except:
		print "Could not run the cdhit program."
		sys.exit()

def dictionary_count(infile):
	current = None
	counts = {}
	for line in Infile:
			if line.startswith('>'):
				line = line.split("\n")[0]
				current = line.split(">")[1]
				counts[current] = 0
			else:
				counts[current] += 1
	return counts


def dictionary_similarhit(directory, infile):
	file = directory + infile
	with open(file, "r+") as Infile:
		current = None
		counts = {}
		for line in Infile:
			if line.startswith('>'):
				line = line.split("\n")[0]
				current = line.split(">")[1]
				counts[current] = 0
			else:
				counts[current] += 1
	return counts
		

def similarhit(directory, infile, counts):
	file = directory + infile
	with open(file, "r+") as Infile:
		list_singleton_trans = []
		for line in Infile:
			if line.startswith('>'):
				line = line.split("\n")[0]
				querry = line.split(">")[1]
				if counts[querry] < 2:
					next_line = Infile.next()
					record = next_line.split(" ")[1]
					record = record.replace("...", "")
					record = record.replace(">", "")
					list_singleton_trans.append(record)
	return list_singleton_trans


print "\n\t 1 Filtering the size of the transcriptome and skimming data..."

"""select the size of the contigs"""
filter_out_dir = "1_filtered_skim_transcriptome_data"
if not os.path.exists(filter_out_dir):
	os.makedirs(filter_out_dir)
skim_filtered = open(filter_out_dir + "/skim_size_filtered.fasta", 'w')
trans_filtered = open(filter_out_dir + "/transcript_size_filtered.fasta", 'w')
for record in SeqIO.parse(skim_contig, 'fasta'):
	length = len(record.seq)
	if  500 <= length <= 1000:
		SeqIO.write(record, skim_filtered, "fasta")
	if length > 1000:
		n = 1
		start = 0
		stop = 1000
		iterations = length/1000 + 1
		for i in range(iterations):
			ID = record.id + "_" + str(n)
			seq = str(record.seq[start:stop])
			if len(seq) > 500:
				new_record = SeqRecord(Seq(seq), id=ID, description="")
				SeqIO.write(new_record, skim_filtered, "fasta")
				stop = stop + 1000
				start = start + 1000
				n = n+1
for record in SeqIO.parse(transciptome, 'fasta'):
	length = len(record.seq)
	if length > 400:
		SeqIO.write(record, trans_filtered, "fasta")
print "\t\tTranscriptome and skimming data size filtered."


print "\n\t 2 Comparing transcriptome and skimming data...that could take a while..."

"""blat between the skimming and transciptome data"""
blat_out_dir = "2_blat_results/"
if not os.path.exists(blat_out_dir):
	os.makedirs(blat_out_dir)
infile_skim = filter_out_dir + "/skim_size_filtered.fasta"
infile_trans = filter_out_dir + "/transcript_size_filtered.fasta"
outfile_blat = blat_out_dir + "skimming_v_transcriptome.pslx"
blat_cmd = "blat -t=dna -q=dna -minIdentity=90 -noHead -out=pslx " + infile_skim + " " + infile_trans + " " + outfile_blat
try:
	with open((blat_out_dir + "blat_screen_out.txt"), 'w') as blat_screen:
		p = subprocess.Popen(blat_cmd,	stdout=blat_screen, shell=True)
		p.communicate()
	print ("\t\tblat completed.")
except:
	print "Could not run the blat program."
	sys.exit()
"""filtering results"""
df = pd.read_csv(outfile_blat, sep="\t", header=None)
df = df.drop_duplicates(subset=9, keep=False)
df_sub = df[16]-df[15] >= df[14]*exon
substantial_hits = list(df[df_sub][9])
print "\t\t" + str(len(substantial_hits)) + " unique transcripts have been selected"
trans_selected = open(blat_out_dir + "trans_selected.fasta", 'w')
for record in SeqIO.parse(transciptome, 'fasta'):
	if record.id in substantial_hits:
		SeqIO.write(record, trans_selected, "fasta")
trans_selected.close()

print "\n\t 3 Remove transcripts with 90% or greater similarity..."

cdhit_out_dir = "3_cd-hit_results/"
if not os.path.exists(cdhit_out_dir):
	os.makedirs(cdhit_out_dir)
infile_cdhit_100_trans = blat_out_dir + "trans_selected.fasta"
outfile_cdhit_100_trans = cdhit_out_dir + "transcripts_hits_cluster_100.fasta"
cdhit(infile_cdhit_100_trans, outfile_cdhit_100_trans, cdhit_out_dir, 1.0, "trans")
outfile_cdhit_90 = cdhit_out_dir + "transcripts_hits_cluster_90.fasta"
cdhit(outfile_cdhit_100_trans, outfile_cdhit_90, cdhit_out_dir, 0.9, "trans")
infile = "transcripts_hits_cluster_90.fasta.clstr"
result = similarhit(cdhit_out_dir, infile, dictionary_similarhit(cdhit_out_dir, infile))
df1 = df[df_sub]
df_mask = df1[df1[9].isin(result)]
df_mask = df_mask.drop_duplicates(subset=13, keep=False)
list_skim = list(df_mask[13])


print "\n\t 4 Remove skim data with 90% or greater similarity..."

skim_selected = open(cdhit_out_dir + "skim_selected.fasta", 'w')
for record in SeqIO.parse(infile_skim, 'fasta'):
	if record.id in list_skim:
		SeqIO.write(record, skim_selected, "fasta")
skim_selected.close()
infile_cdhit_100_skim = cdhit_out_dir + "skim_selected.fasta"
outfile_cdhit_100_skim = cdhit_out_dir + "skim_hits_cluster_100.fasta"
cdhit(infile_cdhit_100_skim, outfile_cdhit_100_skim, cdhit_out_dir, 1.0, "skim")
outfile_cdhit_90_skim = cdhit_out_dir + "skim_hits_cluster_90.fasta"
cdhit(outfile_cdhit_100_skim, outfile_cdhit_90_skim, cdhit_out_dir, 0.9, "skim")
infile = "skim_hits_cluster_90.fasta.clstr"
result = similarhit(cdhit_out_dir, infile, dictionary_similarhit(cdhit_out_dir, infile))

print "\n\t 5 low-copy nuclear markers"

print "\t\t" + str(len(result)) + " low-copy markers have been selected."
df_mask_skim = df_mask[df_mask[13].isin(result)]
list_skim = list(df_mask_skim[13])

markers_out_dir = "4_low_copy_markers/"
if not os.path.exists(markers_out_dir):
	os.makedirs(markers_out_dir)
sequences_for_probe = open((markers_out_dir + "sequences_for_probe.fasta"), 'w')
total_length = 0

for record in SeqIO.parse((cdhit_out_dir + "skim_selected.fasta"), 'fasta'):
	if record.id in list_skim:
		total_length = total_length + len(record.seq)

		SeqIO.write(record, sequences_for_probe, "fasta")
sequences_for_probe.close()
average_length = total_length/len(result)

print "\t\t" + str(total_length) + "bp length of the concatenated low-copy markers."
print "\t\t" + str(average_length) + "bp average length of the low-copy markers."

df1 = (df_mask_skim[16]-df_mask_skim[15])/df_mask_skim[14] 
df2 = df_mask_skim[13]
frames = [df2, df1]
df3 = pd.concat(frames, axis=1)
df3.to_csv((markers_out_dir + "proportion_of_exon_per_markers.csv"), index=False, encoding='utf-8')
df4 = df3.sort_values(0, ascending=True)[[0]]
df4.reset_index(inplace=True)
plt.scatter(df4.index.values, df4[0])
plt.savefig((markers_out_dir + "distribution_exon.png"))


