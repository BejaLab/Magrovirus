import numpy as np
import glob
import argparse
import sys
import time

def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument('dataset', type=str, help="parsed dataset files path")
	parser.add_argument("--outfile", help="Where to write the counts output. (default: stdout)", nargs='?', type=argparse.FileType('w'), default=sys.stdout)
	parser.add_argument("--min_id", help="Cutoff percent-ID. (default: 90)", type=int, default=90)
	return parser.parse_args()

def main():
	start_time = time.time()

	min_id = int(args.min_id)

	dataset = args.dataset

	'''The arrays store: top id and top alignment score for each alignment'''
	ids_array = np.zeros((1, 600000000), dtype=np.int8)
	score_array = np.zeros((1, 600000000), dtype=np.int16)

	# insert in arrays the highest id and score per aln or read
	c = 0
	for file in glob.iglob(dataset + '/*'):
		genome_name = file.split('/')[-1].split('.')[0].split('_')[-1]
		c += 1

		print '\t'.join(map(str,[c, 'P', genome_name, (time.time() - start_time)/60]))
		with open(file,'r') as fh:

			read = None
			pct_id = None
			aln_score = None
			update = 0
			stored_id = None
			stored_score = None
			read_index2 = None

			for line in fh:
				sline = line.split('\t')
				pct_id = int(sline[1])
				if pct_id < min_id:
					continue
				read = int(sline[0])
				aln_score = int(sline[2])

				read_index2 = read - 1

				stored_id = ids_array[0][read_index2]
				stored_score = score_array[0][read_index2]

				if stored_id:# read has been parsed
					if pct_id > stored_id:
						update = 1# update top values
					elif pct_id == stored_id:
						if aln_score > stored_score:# update top values
							update = 1
				elif stored_id == 0:# read has not been parsed
					update = 1
				else:
					continue
				if update:
					ids_array[0][read_index2] = pct_id
					score_array[0][read_index2] = aln_score

	print '\n', c, 'genomes parsed in', (time.time() - start_time)/60, '\n'

	num_genomes = c# parsed genomes
	genome_array = np.zeros((1, num_genomes), dtype=np.uint32)# in this array index_i -> genome_i
	genomes_list = list()

	for file in glob.iglob(dataset + '/*'):
		genome_name = file.split('/')[-1].split('.')[0].split('_')[-1]
		genomes_list.append(genome_name)
		c_genome = len(genomes_list)# counts how many genomes have been processed
		genome_index = c_genome - 1

		print '\t'.join(map(str,[c_genome, 'C', genome_name, (time.time() - start_time)/60]))

		with open(file,'r') as fh:
			read = None
			pct_id = None
			aln_score = None
			update = 0
			stored_id = None
			stored_score = None
			read_index2 = None

			for line in fh:
				sline = line.split('\t')
				pct_id = int(sline[1])
				if pct_id < min_id:
					continue
				read = int(sline[0])
				aln_score = int(sline[2])

				read_index2 = read - 1

				stored_id = ids_array[0][read_index2]
				stored_score = score_array[0][read_index2]

				if (pct_id == stored_id) and (aln_score == stored_score):# if read has the same AS and ID as the top in the array
					genome_array[0][genome_index] += 1
				else:
					continue
	with args.outfile as output_filehandle:
		genome_index = 0# start from the first element stored in genomes_list
		for genome in genomes_list:
			genome_count = genome_array[0][genome_index]
			output_filehandle.write('\t'.join(map(str,[genome, genome_count])) + '\n')
			genome_index += 1

if __name__ == "__main__":
	args = parse_args()
	main()
