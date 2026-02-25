#!/usr/bin/env python

import gzip
import itertools 
import collections
#import urllib 

from pathlib import Path
from urllib import request


#################################################

def load_chrom_namemap_remote(gencode_genome_version):
	ucsc_chrom_to_gencode_chrom = {}
	chrom_namemap_url = "https://raw.githubusercontent.com/dpryan79/ChromosomeMappings/refs/heads/master/{}_UCSC2gencode.txt".format(gencode_genome_version)
	with request.urlopen(chrom_namemap_url) as chrom_namemap_url_file:
		for rawline in chrom_namemap_url_file:
			fields = rawline.decode().strip().split()
			if len(fields)==2:
				ucsc_chrom = fields[0]
				gencode_chrom = fields[1]
				ucsc_chrom_to_gencode_chrom[ucsc_chrom] = gencode_chrom

	return ucsc_chrom_to_gencode_chrom


def load_chrom_namemap_local(chrom_namemap_filename):
	ucsc_chrom_to_gencode_chrom = {}
	with open(chrom_namemap_filename,"rt") as chrom_namemap_file:
		for rawline in itertools.islice(chrom_namemap_file,1,None):
			fields = rawline.strip().split()
			ucsc_chrom = fields[0]
			gencode_chrom = fields[1]
			ucsc_chrom_to_gencode_chrom[ucsc_chrom] = gencode_chrom

	return ucsc_chrom_to_gencode_chrom

#################################################

def load_chrominfo_remote(genome_version):
	chrom_to_size = {}
	chrominfo_url = "https://hgdownload.cse.ucsc.edu/goldenpath/{}/database/chromInfo.txt.gz".format(genome_version)
	with gzip.GzipFile(fileobj=request.urlopen(chrominfo_url)) as chrominfo_url_file:
		for rawbyte in chrominfo_url_file:
			fields = rawbyte.decode().strip().split("\t")
			chrom = fields[0]
			size = int(fields[1])
			chrom_to_size[chrom] = size

	return chrom_to_size


def load_chrominfo_local(chrominfo_filename):
	chrom_to_size = {}
	with open(chrominfo_filename,"rt") as chrominfo_file:
		for rawline in chrominfo_file:
			fields = rawline.strip().split()
			chrom = fields[0]
			size = int(fields[1])
			chrom_to_size[chrom] = size

	return chrom_to_size

#################################################

def load_chrom_fasta_chain(genome_version,ucsc_chrom_to_gencode_chrom,chrom_of_interest_list):
	chrom_fasta_gen_chain = map(lambda ucsc_chrom: load_chrom_fasta(genome_version,ucsc_chrom,ucsc_chrom_to_gencode_chrom[ucsc_chrom]),chrom_of_interest_list)
	for genome_fasta in itertools.chain.from_iterable(chrom_fasta_gen_chain):
		yield genome_fasta 


def load_chrom_fasta(genome_version,ucsc_chrom,gencode_chrom):
	temp_dir = Path("../tmp/") / genome_version
	chrom_fasta_gz_filename = temp_dir / (ucsc_chrom + ".fa.gz")
	with gzip.open(chrom_fasta_gz_filename,"rt") as chrom_fasta_gz_file:
		chrom_fasta = ">" + gencode_chrom + "\n"
		yield chrom_fasta
		for rawline	in itertools.islice(chrom_fasta_gz_file,1,None):
			chrom_fasta = rawline
			yield chrom_fasta

#################################################

def load_blacklist_from_url(genome_version):
	chrom_to_black_list = collections.defaultdict(list)
	blacklist_anno_url = "https://github.com/Boyle-Lab/Blacklist/raw/master/lists/{}-blacklist.v2.bed.gz".format(genome_version)
	with gzip.GzipFile(fileobj=request.urlopen(blacklist_anno_url)) as blacklist_anno_url_file:
		for rawbyte in blacklist_anno_url_file:
			fields = rawbyte.decode().strip().split("\t")
			chrom = fields[0]
			start = int(fields[1])
			end = int(fields[2])
			black_type = fields[3]
			black = chrom,start,end,black_type
			chrom_to_black_list[chrom].append(black)

	return chrom_to_black_list

def load_gap_from_url(genome_version):
	chrom_to_gap_list = collections.defaultdict(list)
	gap_anno_url = "https://hgdownload.cse.ucsc.edu/goldenpath/{}/database/gap.txt.gz".format(genome_version)

	with gzip.GzipFile(fileobj=request.urlopen(gap_anno_url)) as gap_anno_url_file:
		for rawbyte in gap_anno_url_file:
			fields = rawbyte.decode().strip().split("\t")
			chrom = fields[1]
			start = int(fields[2])
			end = int(fields[3])
			gap_type = fields[7]
			gap = chrom,start,end,gap_type
			chrom_to_gap_list[chrom].append(gap)

	return chrom_to_gap_list

#################################################

def load_gencode_gtf_from_url(gencode_version,genome_version):
	genome_version_to_species = {"hg38":"human","mm10":"mouse"}
	species = genome_version_to_species[genome_version]
	gencode_gtf_url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_{}/release_{}/gencode.v{}.primary_assembly.annotation.gtf.gz".format(species,gencode_version,gencode_version)
	#gencode_gtf_url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_{}/release_{}/gencode.v{}.basic.annotation.gtf.gz".format(species,gencode_version,gencode_version)
	with gzip.GzipFile(fileobj=request.urlopen(gencode_gtf_url)) as gencode_gtf_url_file:
		for rawbyte in gencode_gtf_url_file:
			gencode_line = rawbyte.decode().strip()
			yield gencode_line

