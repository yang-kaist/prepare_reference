#!/usr/bin/env python

import gzip

from pathlib import Path

#################################################

def write_chrom_namemap(gencode_chrom_to_ucsc_chrom,sorted_gencode_chrom_list,chrom_namemap_tsv_filename):
	create_parent_directory(chrom_namemap_tsv_filename)
	with open(chrom_namemap_tsv_filename,"wt") as chrom_namemap_tsv_file:
		chrom_namemap_tsv_file.write("UCSC_chrom\tgencode_chrom\n")
		for gencode_chrom in sorted_gencode_chrom_list:
			ucsc_chrom = gencode_chrom_to_ucsc_chrom[gencode_chrom]
			chrom_namemap_tsv_file.write("\t".join([ucsc_chrom,gencode_chrom]) + "\n")

#################################################

def write_chrominfo(chrom_to_size,sorted_chrom_list,chrominfo_filename):
	create_parent_directory(chrominfo_filename)
	with open(chrominfo_filename,"wt") as chrominfo_file:
		for chrom in sorted_chrom_list:
			size = chrom_to_size[chrom]
			chrominfo_file.write( "\t".join(map(str,[chrom,size])) + "\n")

#################################################

def write_bed(bed_interval_list,bed_filename):
	create_parent_directory(bed_filename)
	with open(bed_filename,"wt") as bed_file:
		for bed_interval in bed_interval_list:
			bed_file.write( "\t".join(map(str,bed_interval)) + "\n")

#################################################

def write_gtfgz(gtf_line_list,gtfgz_filename):
	create_parent_directory(gtfgz_filename)
	with gzip.open(gtfgz_filename,"wt") as gtfgz_file:
		for gtf_line in gtf_line_list:
			gtfgz_file.write( gtf_line + "\n")


def write_gtf(gtf_line_list,gtf_filename):
	create_parent_directory(gtf_filename)
	with open(gtf_filename,"wt") as gtf_file:
		for gtf_line in gtf_line_list:
			gtf_file.write( gtf_line + "\n")

#################################################

def create_parent_directory(filename):
	Path.mkdir(Path(filename).parent,parents=True,exist_ok=True)

#################################################