#!/usr/bin/env python

import sys

from genome_toolkit.io import reader
from genome_toolkit.io import writer
from genome_toolkit.io import misc


def main():
	if len(sys.argv) < 3:
		print ("Usage : " + sys.argv[0] + " [chrominfo_filename] [chromsize_bed_filename]")
		sys.exit()
	else:
		chrominfo_filename = sys.argv[1]
		chromsize_bed_filename = sys.argv[2]
		make_chromsize_bed(chrominfo_filename,chromsize_bed_filename)


def make_chromsize_bed(chrominfo_filename,chromsize_bed_filename):
	chrom_to_size = reader.load_chrominfo_local(chrominfo_filename)
	bed_interval_list = misc.chrominfo_to_bed(chrom_to_size)
	writer.write_bed(bed_interval_list,chromsize_bed_filename)


if __name__ == "__main__":
	main()

