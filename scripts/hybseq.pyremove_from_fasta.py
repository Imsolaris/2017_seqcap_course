#!/usr/bin/env python

from Bio import SeqIO
import argparse


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("-i", "--infile", type=str, help="Input file", required=True)
	parser.add_argument("-o", "--out", type=str, help="Outfile", required=True)
	parser.add_argument("-l", "--list_seq", type=str, help="list of sequences to discard", required=True)


	args = parser.parse_args()
	original_file = args.infile
	newfile_name = args.out
	list_sequence = args.list_seq


file = open(original_file, "r")
newfile = open(newfile_name, 'w')
with open(list_sequence) as f:
        remove_list = f.read().splitlines()
for record in SeqIO.parse(file, 'fasta'):
        if record.id not in remove_list:
                SeqIO.write(record, newfile, "fasta")
file.close()
