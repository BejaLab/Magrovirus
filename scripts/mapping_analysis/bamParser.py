"""
requires HTseq (http://www-huber.embl.de/users/anders/HTSeq/)
parse SAM/BAM alignment and extract alignment hits by percent-identity and alignment score.
requires SAM/BAM file to have an MD string
and (http://code.google.com/p/pysam/).
Input is a bam/sam file sorted or a stream.
"""

import os, sys
import re
import HTSeq
import itertools
import argparse
import numpy as np

current_dir = os.getcwd()

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('input_file', type=str, help="Required: <input file> , for stdin use '-'")
	parser.add_argument("-i", "--min_id", help = "minimal id percent of hit to be considered. (default: 80)", default=80, type=int)
	parser.add_argument("-l", "--min_len", help = "minimal read length to be considered. (default: 80)", default=80, type=int)
	parser.add_argument("-c", "--max_clip", help = "proportion of bases clipped from read for alignment. (default: 0.3)", default=0.3, type=float)
	parser.add_argument("--outfile", help="Where to write the counts output. (default: stdout)", nargs='?', type=argparse.FileType('w'), default=sys.stdout)
	parser.add_argument("--bam_file", help="Path to the bam file, required to catch false empty BAM files errors", type=str, required=True)

	args = parser.parse_args()

	min_read_len  = args.min_len
	min_id = args.min_id
	max_clip_ = args.max_clip

	paired_end = True
	pe_order = 'n'

	if args.input_file:
		try:
			if args.input_file == '-':
				seqfile = HTSeq.SAM_Reader(sys.stdin)
				try:
					if paired_end:
						read_seq_iter = iter(seqfile)
						first_read = read_seq_iter.next()
						if pe_order == 'p':
							reader = HTSeq.pair_SAM_alignments_with_buffer(seqfile)
						elif pe_order == 'n':
							reader = HTSeq.pair_SAM_alignments(seqfile)
					else:
						reader = seqfile
				except StopIteration as e:
					if os.path.os.path.getsize(args.bam_file) < 400:
						error_message = "Empty SAM/BAM file.\n"
						sys.stderr.write(error_message)
						sys.stderr.write(str(e)+'\n')
						sys.exit(0)
					else:
						error_message = "failed processing SAM/BAM file. Try submission again."
						sys.stderr.write(str(e)+'\n')
						sys.exit(error_message)
					sys.stderr.write(str(e)+'\n')
					raise
			elif args.input_file != '-':
					seqfile= HTSeq.SAM_Reader(args.input_file)
					if paired_end:
						read_seq_iter = iter(seqfile)
						first_read = read_seq_iter.next()
						read_seq = itertools.chain( [ first_read ], read_seq_iter)
						reader = HTSeq.pair_SAM_alignments(read_seq)
						if pe_order == 'p':
							reader = HTSeq.pair_SAM_alignments_with_buffer(reader)
						elif pe_order == 'n':
							reader = HTSeq.pair_SAM_alignments(reader)
					else:
						reader = seqfile
		except Exception as e:
			error_message = "failed processing SAM/BAM file."
			sys.stderr.write(str(e))
			sys.exit(error_message)
			raise

	read_counter = 0
	dataset_id_dict = dict()

	ids_array = np.zeros((1, 300000000),dtype=np.int8)# stores the highest id
	score_array = np.zeros((1, 300000000),dtype=np.int16)# stores the sum of the bowtie2 alignment score
	genome = ''

	for alignment in reader:
		read_counter += 1
		iv_seq_good_1 = False
		iv_seq_good_2 = False
		aln_score = 0

		if  genome == '':
			try:
				for aln_i in alignment:
					genome = aln_i.iv.chrom
					continue
			except:
				pass
		try:
			read_1_name, percent_1_id, iv_seq_good_1, aln1_score = calc(alignment[0], min_id, max_clip_, min_read_len)
			read_index = read_1_name - 1
		except:
			iv_seq_good_1 = False
			aln1_score = 0

			pass
		try:
			read_2_name, percent_2_id, iv_seq_good_2, aln2_score = calc(alignment[1], min_id, max_clip_, min_read_len)
			read_index = read_2_name - 1
		except:
			iv_seq_good_2 = False
			aln2_score = 0

			pass

		aln_score = aln1_score + aln2_score

		if iv_seq_good_1:
			aln_id = percent_1_id
			if iv_seq_good_2:
				if percent_1_id < percent_2_id:
					aln_id = percent_2_id
		elif iv_seq_good_2:
			aln_id = percent_2_id
		else:
			continue

		if ids_array[0][read_index] == 0:
			ids_array[0][read_index] = aln_id
			score_array[0][read_index] = aln_score
		elif aln_id >= ids_array[0][read_index]:
			if aln_score >= score_array[0][read_index]:
				ids_array[0][read_index] = aln_id
				score_array[0][read_index] = aln_score

	with args.outfile as outputFileHandle:
		for readindex in np.nonzero(ids_array)[1]:
			outputFileHandle.write('\t'.join(map(str,[readindex + 1, ids_array[0][readindex], score_array[0][readindex]]))+'\n')

def calc(aln, min_id, max_clip_, min_read_len):
	if (aln is not None) and (aln.aligned):
		read_name = aln.read.name
		read_name = int(read_name.split('.')[1])
		read_length = len(aln.read.seq)
		opt_fields = aln.optional_fields
		if min_id > 0:
			cigar_string = parse_cigar(aln.original_sam_line.split('\t')[5])
			cigar_soft_clipped, cigar_M, cigar_insertions, cigar_deletions, cigar_insertions = parse_cigar_alignment(cigar_string)
			score, md_matches, md_deletions, md_mismatches = parse_opt_fields(opt_fields)
			clipped = (float(cigar_soft_clipped)/float(read_length))
			percent_id = int((100.0 * ((float(md_matches) / (float(read_length - cigar_soft_clipped + cigar_insertions + cigar_deletions))))))
			if percent_id >= min_id:
				if (float(cigar_soft_clipped)/float(read_length)) <= float(max_clip_):
					iv_seq_good = True
					id_check = True
				else:
					iv_seq_good = False
					return read_name, percent_id, iv_seq_good, score
			else:
				id_check = False
				iv_seq_good = False
				return read_name, percent_id, iv_seq_good, score
		else:
			id_check = True
			iv_seq_good = True
		if read_length >= min_read_len:
			if id_check == True:
				iv_seq_good = True
			else:
				iv_seq_good = False
				return read_name, percent_id, iv_seq_good, score
		else:
			iv_seq_good = False
			return read_name, percent_id, iv_seq_good, score
	else:
		iv_seq_good = False
		return '', 0, iv_seq_good, 0

	return read_name, percent_id, iv_seq_good, score

def parse_cigar_alignment(cigar_string):
	# parse cigar string to get alignment data: total number of insertions, clipped bases, deletions, etc
	cigar_soft_clipped = 0
	cigar_M = 0
	cigar_insertions = 0
	cigar_deletions = 0
	cigar_insertions = 0
	for c in cigar_string:
		if 'S' in c:
			c = c.replace('S','')
			cigar_soft_clipped += int(c)
		elif 'M' in c:
			c = c.replace('M','')
			cigar_M += int(c)
		elif 'D' in c:
			c = c.replace('D','')
			cigar_deletions += int(c)
		elif 'I' in c:
			c = c.replace('I','')
			cigar_insertions += int(c)
	return cigar_soft_clipped, cigar_M, cigar_insertions, cigar_deletions, cigar_insertions

def parse_opt_fields(opt_fields):
	# parse optional fields column in sam file.
	# get MD string alignment data.
	# get score --> AS
	score = int()
	md_matches = int()
	md_deletions = int()
	md_mismatches = int()
	for field in opt_fields:
		if 'MD' in field[0]:
			md_list = re.split('(\D|\W)', field[1])
			for i in md_list:
				if i:
					if i.isdigit():
						md_matches += int(i)
					elif "^" in i:
						md_deletions += 1
					else:
						md_mismatches += 1
			total_md = md_mismatches + md_deletions + md_matches
		if 'AS' in field[0]:
			score = field[1]
	return score, md_matches, md_deletions, md_mismatches

def parse_cigar(s):
	s = s
	l = list()
	ll = list()
	for i in s:
		if i.isdigit():
			ll.append(i)
		else:
			ll.append(i)
			l.append(''.join(ll))
			ll = list()
	return l

main()
