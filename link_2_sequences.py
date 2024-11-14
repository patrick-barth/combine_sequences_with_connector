#!/usr/bin/env python

import argparse
import os
from Bio import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import pprint
import random
from random_sequence_generator import generate_random_seq

##############################################
# Reads 2 fasta files with sequences and     #
#  connects them with either a random        #
#  sequence or a supplied sequence           #
##############################################

__author__ = "Patrick Barth" 
__email__ = "Patrick.Barth@biologie.uni-regensburg.de"

parser = argparse.ArgumentParser()
parser.add_argument('--5_prime_file', '-5', default=False, dest='five_prime_file',
					help='File containing sequences to be put in 5\' end')
parser.add_argument('--3_prime_file', '-3', default=False, dest='three_prime_file',
					help='File containing sequences to be put in 3\' end')
parser.add_argument('--connector_seq', '-c', default=False,
					help='File containing sequences to connect 5\' and 3\' sequences')
parser.add_argument('--illegal_seq_file', '-i', default=False,
					help='File containing motifs not to used. Omits sequences motifs them and does not build them in randomly generated sequences.')
parser.add_argument('--omit_standard_seq', '-s', action=argparse.BooleanOptionalAction, default=False,
					help='Removes standard motifs')
parser.add_argument('--number_seqs','-n', default='10',
					help='Number of sequences to be generated')
parser.add_argument('--use_same_connector','-a', action='store_true', default=False,
                    help='If true only one connector sequence is generated and used. Only applies when a random connector sequence is generated'),
parser.add_argument('--connector_length','-b', default=100,
                    help='Length of connector sequence. Only applies when random connector is generated'),
parser.add_argument('--seq_gc','-g', default='50',
					help='Chance for a base to be G or C. Only applies if random linker is used.')
parser.add_argument('--max_len', '-l', default='100000',
					help='Length of generated sequences. If combined sequences are longer than this they are not generated.'),
parser.add_argument('--enforce_longer_sequences','-e',action='store_true', default=False,
                    help='Enforces the output of sequences longer than the chosen maximum.'),
parser.add_argument('--output','-o', default='./output',
					help='File to write output to')
args = parser.parse_args()

# Print overview of used parameters
print("Used parameters:")
pprint.pprint(vars(args))



###################
### main script ###
###################

def main(five_prime_file,three_prime_file,connector_file,use_same_connector,connector_length,seq_gc,number_seqs,omit_standard_seq,output):
	five_prime_seqs	= import_fasta(five_prime_file)
	three_prime_seqs = import_fasta(three_prime_file)
	connector_length = int(connector_length)
	seq_gc = float(seq_gc)
	number_seqs = int(number_seqs)

	connector_seqs = []
	merged_seqs = []

	if connector_file:
		connector_seqs = import_fasta(connector_file)
	else:
		if use_same_connector:
			connector_seqs = generate_random_seq(connector_length,seq_gc,omit_standard_seq)
		else:
			for i in range(number_seqs):
				connector_seqs.append(generate_random_seq(connector_length,seq_gc,omit_standard_seq))

	for i in range(number_seqs):
		# Get random entry from sequences
		# TODO: Start here if specific sequences are to be picked
		five_prime = random.choice(five_prime_seqs)
		connector = random.choice(connector_seqs)
		three_prime = random.choice(three_prime_seqs)
		
		connector_id = connector.id if connector_file else f"connector_{str(i+1)}"
		connector_seq = connector.seq if connector_file else connector

		combined_seq = five_prime.seq + connector_seq + three_prime.seq

		ID = f"{five_prime.id}_{connector_id}_{three_prime.id}"
		positions = f"pos_5:1-{str(len(five_prime.seq))}"
		positions = positions + f"pos_con:{str(len(five_prime.seq)+1)}-{str(len(five_prime.seq) + len(connector))}"
		positions = positions + f"pos_3:{str(len(five_prime.seq) + len(connector) + 1)}-{len(combined_seq)}"


		
		combined_seq_record = SeqRecord(combined_seq, id=ID, description=positions)
		merged_seqs.append(combined_seq_record)

	with open(output, "w") as output_file:
		SeqIO.write(merged_seqs, output_file, "fasta")



#######################
###    functions    ###
#######################
	
def import_fasta(path):
	collected_seqs = []
	for record in SeqIO.parse(path, "fasta"):
		collected_seqs.append(record)
	return collected_seqs



#########################
### start main script ###
#########################

if __name__== "__main__":
	main(five_prime_file = args.five_prime_file,
		three_prime_file = args.three_prime_file,
		connector_file = args.connector_seq,
		use_same_connector = args.use_same_connector,
		connector_length = args.connector_length,
		seq_gc = args.seq_gc,
		number_seqs = args.number_seqs,
		omit_standard_seq = args.omit_standard_seq,
		output = args.output)