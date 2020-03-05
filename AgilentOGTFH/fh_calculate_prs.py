"""
A script to calculate a Polygenic Risk Score (PRS) for the FH assay

See - https://www.ncbi.nlm.nih.gov/pubmed/23433573 (talmud 2013) for info on how calculation works.

Usage:

python fh_calculate_prs.py \
--genotypes "$seqId"_"$sampleId"_snps_fixed_decomposed.csv \
--annotations /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_snp_annotations.yaml \
--output_name "$seqId"_"$sampleId" \
--sample_id "$sampleId"

Input:
--genotypes - csv file containing the genotypes of the snps
--annotations - yaml config file with weights etc
--output_name - what to call the output
--sample_id - sample id

Output:

*_snps.csv - csv file containing all the genotyped snps and their annotations
*_snps_prs.csv - csv file containing just the snps used for the PRS
*_prs_score.csv - csv file containing the PRS score and which decile the sample falls in.


"""

import argparse
import pandas as pd
import yaml
import csv



# parse args
parser = argparse.ArgumentParser(description='')
parser.add_argument('--genotypes', type=str, nargs=1, required=True,
				help='Input CSV with genotypes')
parser.add_argument('--annotations', type=str, nargs=1, required=True,
				help='Input YAML with annotations for the SNPs')
parser.add_argument('--output_name', type=str, nargs=1, required=True,
				help='Where to put the output')
parser.add_argument('--sample_id', type=str, nargs=1, required=True,
				help='The sample ID')

args = parser.parse_args()


# utility functions

def parse_config(config_location):
	"""
	Parse the YAML config file.
	"""
	with open(config_location, 'r') as stream:
		return yaml.safe_load(stream)
	
def annotate_genotype_df(df, snps_dict, key):
	"""
	Function to get annotation from a dictionary \
	and add to genotype df.

	"""
	
	chrom = df['CHROM']
	
	pos = df['POS']
	
	snp_key = f'{chrom}:{pos}'
		
	annotation = snps_dict[snp_key][key]
	
	return annotation

def count_risk_alleles(df):
	"""
	Count number of risk alleles in a genotype.

	"""
	
	gt = df[f'{sample_id}.GT']
	risk_allele = df['risk_allele']
	
	if '/' in gt:
		
		gt = gt.split('/')
		
	elif '|' in gt:
		
		gt = gt.split('|')
		
	else:
		
		raise Exception('oh no')
		
	return gt.count(risk_allele)

def get_apoe_snp_score(snps_for_prs_apoe, sample_id):
	"""
	Calculate the PRS using the APOE SNPs.

	"""

	print(snps_for_prs_apoe[snps_for_prs_apoe['snp_id'] =='rs429358'][f'{sample_id}.GT'])
	
	rs429358_gt = snps_for_prs_apoe[snps_for_prs_apoe['snp_id'] =='rs429358'][f'{sample_id}.GT'].iloc[0]
	
	if '/' in rs429358_gt:
		
		rs429358_gt = rs429358_gt.split('/')
		
	elif '|' in rs429358_gt:
		
		rs429358_gt = rs429358_gt.split('|')
		
	else:
		
		raise Exception('oh no')
		
	for gt in rs429358_gt:
		
		if gt not in ['C', 'T']:
			
			raise Exception('Invalid genotype for rs429358')
			
	
	rs429358_t_count = rs429358_gt.count('T')
	
	rs7412_gt = snps_for_prs_apoe[snps_for_prs_apoe['snp_id'] =='rs7412'][f'{sample_id}.GT'].iloc[0]
	
	if '/' in rs7412_gt:
		
		rs7412_gt = rs7412_gt.split('/')
		
	elif '|' in rs7412_gt:
		
		rs7412_gt = rs7412_gt.split('|')
		
	else:
		
		raise Exception('An invalid genotype')
		
	for gt in rs7412_gt:
		
		if gt not in ['C', 'T']:
			
			raise Exception('Invalid genotype for rs7412')
	
	rs7412_t_count = rs7412_gt.count('T')
	
	if rs429358_t_count == 2 and rs7412_t_count == 2:
		
		return -0.9
	
	elif rs429358_t_count == 2 and rs7412_t_count == 1:
		
		return -0.4
	
	elif rs429358_t_count == 1 and rs7412_t_count == 1:
		
		return -0.2
	
	elif rs429358_t_count == 2 and rs7412_t_count == 0: 

		return 0.0
	
	elif rs429358_t_count == 1 and rs7412_t_count == 0:  
	
		return 0.1
	
	elif rs429358_t_count == 0 and rs7412_t_count == 0: 
		
		return 0.2
	
	else:
		
		raise Exception('An invalid genotype combination for the APOE genes.')
		
def get_decile(score):
	"""
	Workout which decile the score falls in.

	ADD PAPER REF

	"""
	
	if score >= -0.5 and score < 0.58:
		
		return 1, 'Low'
	
	elif score >= 0.58 and score < 0.73:
		
		return 2, 'Low'
	
	elif score >= 0.73 and score < 0.81:
		
		return 3, 'Low'
	
	elif score >= 0.81 and score < 0.88:
		
		return 4, 'Borderline'
	
	elif score >= 0.88 and score < 0.93:
		
		return 5, 'Borderline'
	
	elif score >= 0.93 and score < 0.98:
		
		return 6, 'High'
	
	elif score >= 0.98 and score < 1.02:
		
		return 7, 'High'
	
	elif score >= 1.02 and score < 1.08:
		
		return 8, 'High'
	
	elif score >= 1.08 and score < 1.16:
		
		return 9, 'High'
	
	elif score >= 1.16 and score < 1.46:
		
		return 10, 'High'
	
	else:
		
		return 'NA', 'NA'


def check_qc(df, min_depth, min_gq, sample_id):
	"""
	Check that we have high enough depth and gq

	"""
	
	qc = []
	
	depth = df[f'{sample_id}.NR']
	
	gq = df[f'{sample_id}.GQ']
	
	if depth < min_depth:
		
		qc.append('low_depth')
		
	if depth < min_gq:
		
		qc.append('low_gq')
		
	if len(qc) == 0:
		
		return 'pass'
	
	else:
		
		return '|'.join(qc)


# Main Program

sample_id = args.sample_id[0]
min_gq = 90
min_depth = 20

# TSV file containing the genotypes for the samples.
genotype_df = pd.read_csv(args.genotypes[0], sep='\t')

# YAML file containing the annotations such as rs number, weights etc.
snp_info_dict = parse_config(args.annotations[0])

if genotype_df.shape[0] < 20:

	print ('Could not read CSV (expecting more rows)')

else:

	# Add annotations
	genotype_df['snp_id'] = genotype_df.apply(annotate_genotype_df, axis=1, args=(snp_info_dict, 'id',))
	genotype_df['reason'] = genotype_df.apply(annotate_genotype_df, axis=1, args=(snp_info_dict, 'reason',))
	genotype_df['gene'] = genotype_df.apply(annotate_genotype_df, axis=1, args=(snp_info_dict, 'gene',))
	genotype_df['is_prs'] = genotype_df.apply(annotate_genotype_df, axis=1, args=(snp_info_dict, 'is_prs',))
	genotype_df['is_apoe_snp'] = genotype_df.apply(annotate_genotype_df, axis=1, args=(snp_info_dict, 'is_apoe_snp',))
	genotype_df['reference_reads'] = genotype_df[f'{sample_id}.NR'] - genotype_df[f'{sample_id}.NV']
	genotype_df['alt_reads'] = genotype_df[f'{sample_id}.NV']
	genotype_df['qc'] = genotype_df.apply(check_qc, axis=1,args=(min_depth, min_gq, sample_id))

	# Save all SNPs to CSV for viewing by user
	genotype_df[['CHROM',
	 'POS',
	 'REF',
	 'ALT',
	 'snp_id',
	 'gene',
	  f'{sample_id}.GT',
	  f'{sample_id}.GQ',
	  'qc',
	  'reference_reads',
	  'alt_reads',
	  'reason']].to_csv(args.output_name[0] + '_snps.csv', sep='\t', index=False)


	# get snps for PRS but exclude those in APOE as they are treated differently
	snps_for_prs_non_apoe = genotype_df[(genotype_df['is_prs']==True) & (genotype_df['is_apoe_snp']==False)]

	# Annotate the prs snps
	snps_for_prs_non_apoe['risk_allele'] = snps_for_prs_non_apoe.apply(annotate_genotype_df, axis=1, args=(snp_info_dict, 'risk_allele'))
	snps_for_prs_non_apoe['prs_weight'] = snps_for_prs_non_apoe.apply(annotate_genotype_df, axis=1, args=(snp_info_dict, 'weight'))


	# Add count of risk alleles and the count of risk alleles multiplied by the weight
	snps_for_prs_non_apoe['risk_allele_count'] = snps_for_prs_non_apoe.apply(count_risk_alleles, axis=1)
	snps_for_prs_non_apoe['risk_allele_count'] = snps_for_prs_non_apoe['risk_allele_count'].astype(int)
	snps_for_prs_non_apoe['per_snp_prs_score'] = snps_for_prs_non_apoe['prs_weight'] * snps_for_prs_non_apoe['risk_allele_count']

	# add up all the scores
	non_apoe_score = snps_for_prs_non_apoe['per_snp_prs_score'].sum()

	# now do the APOE SNPs
	snps_for_prs_apoe = genotype_df[genotype_df['is_apoe_snp']==True]
	apoe_score = get_apoe_snp_score(snps_for_prs_apoe, sample_id)
	final_score = apoe_score + non_apoe_score

	# merge the two dataframes for saving to disk for the user

	merged_df = snps_for_prs_non_apoe.append(snps_for_prs_apoe, sort=True)

	# save to disk
	merged_df[['CHROM',
	 'POS',
	 'REF',
	 'ALT',
	 'snp_id',
	 'gene',
	 f'{sample_id}.GT',
	 'risk_allele',
	 'reference_reads',
	 'alt_reads',
	 'qc',
	 'prs_weight',
	 'risk_allele_count',
	 'per_snp_prs_score' ]].to_csv(args.output_name[0] + '_snps_prs.csv', sep='\t', index=False)

	# calculate decile

	decile, decile_description = get_decile(final_score)


	# check all snps in merged are pass qc status

	unique_qcs = merged_df['qc'].unique()

	prs_snps_qc = 'unknown'

	if unique_qcs[0] == 'pass' and len(unique_qcs) ==1:

		prs_snps_qc = 'pass'

	else:

		prs_snps_qc = 'fail'

	with open(args.output_name[0] + '_prs_score.csv', mode='w') as prs_score:

		prs_writer = csv.writer(prs_score, delimiter='\t')
		prs_writer.writerow(['score', 'decile', 'decile_description', 'qc'])
		prs_writer.writerow([round(final_score,5), decile, decile_description, prs_snps_qc ])