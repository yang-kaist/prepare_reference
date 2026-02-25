#!/usr/bin/env python

import itertools 
import collections 

from pathlib import Path
from urllib import request

#################################################

def download_chrom_fasta_gz(genome_version,chrom_of_interest_list):
	temp_dir = Path("../tmp/") / genome_version
	Path.mkdir(temp_dir, parents=True, exist_ok=True)
	for chrom in chrom_of_interest_list:
		chrom_sequence_file_url = "https://hgdownload.soe.ucsc.edu/goldenPath/{}/chromosomes/{}.fa.gz".format(genome_version,chrom)
		chrom_fasta_gz_filename = temp_dir / (chrom + ".fa.gz")
		request.urlretrieve(chrom_sequence_file_url,chrom_fasta_gz_filename)

def clear_tmp(genome_version,chrom_of_interest_list):
	temp_dir = Path("../tmp/") / genome_version
	for temp_filename in temp_dir.glob("*"):
		Path(temp_filename).unlink(missing_ok=True)

#def clear_tmp(genome_version,chrom_of_interest_list):
#	temp_dir = Path("../tmp/") / genome_version
#	for chrom in chrom_of_interest_list:
#		chrom_fasta_gz_filename = temp_dir / (chrom + ".fa.gz")
#		Path(chrom_fasta_gz_filename).unlink(missing_ok=True)

#################################################

def download_gencode_gtf_gz(genome_version,gencode_version):
	temp_dir = Path("../tmp/") / genome_version
	Path.mkdir(temp_dir, parents=True, exist_ok=True)
	genome_version_to_species = {"hg38":"human","hg19":"human","mm10":"mouse","mm39":"mouse"}
	species = genome_version_to_species[genome_version]
	#gencode_gtf_url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_{}/release_{}/gencode.v{}.primary_assembly.annotation.gtf.gz".format(species,gencode_version,gencode_version)
	gencode_gtf_url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_{}/release_{}/gencode.v{}.basic.annotation.gtf.gz".format(species,gencode_version,gencode_version)
	gencode_gtf_gz_filename = temp_dir / "gencode.gtf.gz"
	request.urlretrieve(gencode_gtf_url,gencode_gtf_gz_filename)

#################################################

def sort_chrom_list(chrom_list):
	#predefined_sex_chrom_slim_set = set(["X"])
	predefined_sex_chrom_slim_set = set(["X","Y"])
	predefined_mito_chrom_slim_set = set(["M","MT","Mt"])

	autosome_slim_set = set()
	sex_chrom_slim_set = set()
	mito_chrom_slim_set = set()
	other_chrom_set = set()
	for chrom in chrom_list:
		chrom_slim = chrom.replace("chr","")
		if chrom_slim.isnumeric():
			autosome_slim_set.add(int(chrom_slim))
		elif chrom_slim in predefined_sex_chrom_slim_set:
			sex_chrom_slim_set.add(chrom_slim)
		elif chrom_slim in predefined_mito_chrom_slim_set:
			mito_chrom_slim_set.add("M")
		else:
			other_chrom_set.add(chrom)

	autosome_list = ["chr"+str(chrom_slim) for chrom_slim in sorted(autosome_slim_set) ]
	sex_chrom_list = ["chr"+str(chrom_slim) for chrom_slim in sorted(sex_chrom_slim_set) ]
	mito_chrom_list = ["chr"+str(chrom_slim) for chrom_slim in sorted(mito_chrom_slim_set) ]
	other_chrom_list = sorted(other_chrom_set)
	sorted_chrom_list = autosome_list + sex_chrom_list + mito_chrom_list + other_chrom_list
	return sorted_chrom_list

#################################################

def filter_chrom_list(chrom_list):
	predefined_sex_chrom_slim_set = set(["X"])
	filtered_chrom_list = []
	for chrom in chrom_list:
		chrom_slim = chrom.replace("chr","")
		if chrom_slim.isnumeric() or chrom_slim in predefined_sex_chrom_slim_set:
			filtered_chrom_list.append(chrom)

	return filtered_chrom_list

def is_autox(chrom):
	chrom_slim = chrom.replace("chr","")
	if chrom_slim.isnumeric() or chrom_slim == "X":
		return True
	else:
		return False

#################################################

def chrominfo_to_bed(chrom_to_size):
	interval_list = []
	for chrom,size in chrom_to_size.items():
		interval = (chrom,0,size)
		interval_list.append(interval)

	return interval_list


def chain_interval_list_by_chrom_order(chrom_to_interval_list,chrom_list):
	for chrom in chrom_list:
		if chrom in chrom_to_interval_list:
			interval_list = chrom_to_interval_list[chrom]
			for interval in interval_list:
				yield interval

#################################################

def swap_dict(key_to_value):
	if len(set(key_to_value.keys())) != len(set(key_to_value.values())):
		print("key or value is not unique.")

	value_to_key = {}
	for key in key_to_value.keys():
		value = key_to_value[key]
		value_to_key[value] = key

	return value_to_key

#################################################

def get_commentation_line(rawline_gen,n_max,comment_char):
	for rawline in itertools.islice(rawline_gen,n_max):
		if rawline.startswith(comment_char):
			yield rawline

#################################################
