import pandas as pd
import glob
from pysam import VariantFile
import csv
import argparse

"""
Run in run folder like:

python /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/sv_calling/create_cnv_report.py \
        --gene_bed_file /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI_b37_CNV.bed \
        --bam_list HighCoverageBams.list \
        --worksheet_id $seqId \
        --output sv_analysis/"$seqId"_cnvReport.csv
"""

# Get arguments
parser = argparse.ArgumentParser(
	formatter_class=argparse.RawTextHelpFormatter,
	description='Program to merge SV calls from Manta and Decon.')

parser.add_argument('--gene_bed_file', type=str, nargs=1, required=True,
					help='The bed file of regions with their associated gene annotations.')

parser.add_argument('--bam_list', type=str, nargs=1, required=True,
					help='A text file containing the high coverage bams used by Manta and Decon.')

parser.add_argument('--worksheet_id', type=str, nargs=1, required=True,
					help='The worksheet id.')

parser.add_argument('--output', type=str, nargs=1, required=True,
					help='Name of the output CSV.')

args = parser.parse_args()

# config variables

gene_bed_file = config = args.gene_bed_file[0]
bam_list = args.bam_list[0]
worksheet_id = args.worksheet_id[0]
output = args.output[0]

decon_csv = f'old/{worksheet_id}*_all.txt'
manta_path = 'AgilentOGTFH/*/*_sv_filtered.vcf.gz'
sample_depths = 'AgilentOGTFH/*/*_DepthOfCoverage.sample_summary'


# Useful Functions

def create_manta_df(manta_vcfs, high_cov_sample_list):
	"""
	Loop through manta VCFs and create a dataframe.

	"""
	
	manta_list = []
	
	for vcf in manta_vcfs:
		
		sample_name = vcf.split('/')[-2]
		
		# ignore vcfs made from non high coverage samples.
		if sample_name in high_cov_sample_list:
		
			vcf_in = VariantFile(vcf)

			for rec in vcf_in.fetch():

				chrom = rec.chrom
				start = rec.pos
				end = rec.stop
				ref = rec.ref
				alt = rec.alts[0]
				sv_type = rec.info['SVTYPE']
				sv_filter = ';'.join(rec.filter.keys())

				manta_list.append([sample_name, 'manta', chrom, start, end, ref, alt, sv_type, '-', sv_filter, 'NA', 'NA', 'NA'])
			
			
	manta_df = pd.DataFrame(manta_list, columns=['sample_name',
												 'caller',
												 'chromosome',
												 'start',
												 'end',
												 'ref',
												 'alt',
												 'variant_type',
												 'n_affected_exons',
												 'variant_filter',
												 'decon_correlation',
												 'exon_start',
												 'exon_end'])
	
	return manta_df


def get_overlapping_genes(df, grouped_bed):
	"""
	Get the custom roi file and check to see if any of our SVs overlap with the regions.
	"""

	sv_start = df['start']

	if sv_start == 'NA':

		return 'NA'
	
	if pd.isna(df['gene']) == False:
		
		return df['gene']

	sv_chrom = str(df['chromosome'])
	sv_start = int(df['start'])
	sv_end = int(df['end'])

	gene_list = []

	for gene in grouped_bed.itertuples():

		gene_chrom = str(gene.chrom)
		gene_start = int(gene.start)
		gene_end = int(gene.end)

		# variant start is within the gene
		if gene_chrom == sv_chrom and gene_start <= sv_start <= gene_end:

			gene_list.append(gene.gene)

		# variant end is within the gene
		elif gene_chrom == sv_chrom and gene_start <= sv_end <= gene_end:

			gene_list.append(gene.gene)

		# variant start is before gene and variant end is after gene e.g. variant engulfs gene.
		elif gene_chrom == sv_chrom and (sv_start <= gene_start and sv_end >= gene_end):

			gene_list.append(gene.gene)

	if len(gene_list) > 0:

		return '|'.join(list(set(gene_list)))

	else:

		return 'NA'

def add_sample_filter(df):
	"""
	Add a sample level filter for low coverage samples and \
	samples with low r2 score from Decon.

	"""
	
	sample_filter = []
	
	if df['decon_correlation'] != 'NA':
		
		if df['decon_correlation'] < 0.98:
			
			sample_filter.append('R2<0.98')
		
	if df['mean'] != 'NA':
				
		if df['mean'] < 160:
					
			sample_filter.append('Depth<160')
	
	if len(sample_filter) ==0:
		
		return 'PASS'
	
	else:
		
		return ';'.join(sample_filter)       

def fix_variant_type(df):
	"""
	Unify the variant type field between Manta and Decon.

	"""
	if df['variant_type'] == 'BND':
		
		return 'breakend'
	
	elif df['variant_type'] == 'INS':
		
		return 'insertion'
	
	elif df['variant_type'] == 'DUP':
		
		return 'duplication'
	
	elif df['variant_type'] == 'DEL':
		
		return 'deletion'
	
	elif df['variant_type'] == 'INV':
		
		return 'inversion'
	
	else:
		
		return df['variant_type']


# Open high coverage bam list and create list of high coverage bam files.
high_cov_sample_list = []

with open(bam_list, newline='') as csvfile:
	spamreader = csv.reader(csvfile, delimiter=' ')
	for row in spamreader:
		high_cov_sample_list.append(row[0].split('/')[-2])


decon_csvs = glob.glob(decon_csv)


if len(decon_csvs) == 0:

	master_decon_csv = pd.DataFrame(columns =['sample_name',
		'caller',
		'chromosome',
		'start',
		'end',
		'ref',
		'variant_type',
		'n_affected_exons',
		'variant_filter',
		'gene',
		'decon_correlation',
		'exon_start',
		'exon_end' ] )
	print ('Could not find CNV csv file - decon does not create file if no CNVs are called so ignoring problem.')


# Open decon CSV as dataframe
else:

	master_decon_df = pd.DataFrame()

	for decon_csv in decon_csvs:

		decon_df = pd.read_csv(decon_csv, sep='\t', na_filter = False)
		# Change column names and add other columns that we need
		decon_df['correlation'] = decon_df['Correlation']
		decon_df['chromosome'] = decon_df['Chromosome']
		decon_df['start'] = decon_df['Start']
		decon_df['end'] = decon_df['End']
		decon_df['variant_type'] = decon_df['CNV.type']
		decon_df['variant_filter'] = 'PASS'
		decon_df['exons_affected'] = decon_df['N.exons']
		decon_df['gene'] = decon_df['Gene']
		decon_df['caller'] = 'decon'
		decon_df['sample_name'] = decon_df['Sample'].apply(lambda x: x.split('_')[-1])
		decon_df['decon_correlation'] = decon_df['Correlation']
		decon_df['n_affected_exons'] = decon_df['N.exons']
		decon_df['ref'] = 'NA'
		decon_df['alt'] = 'NA'


		if 'Custom.first' in decon_df.columns:

			decon_df['exon_start'] = decon_df['Custom.first']
			decon_df['exon_end'] = decon_df['Custom.last']

		else:

			decon_df['exon_start'] = 'NA'
			decon_df['exon_end'] = 'NA'

		master_decon_df = master_decon_df.append(decon_df)

	master_decon_df = master_decon_df.drop_duplicates(['sample_name', 'start', 'end', 'variant_type'])

# Get the list of manta vcfs to open and then create Manta df
manta_vcfs = glob.glob(manta_path)
manta_df = create_manta_df(manta_vcfs, high_cov_sample_list)

# Get list of sample depth summary files to open
sample_depths_files = glob.glob(sample_depths)

# Merge the manta and decon dfs
merged_df = manta_df.append(master_decon_df, sort=False)[['sample_name',
 'caller',
 'chromosome',
 'start',
 'end',
 'ref',
 'alt',
 'variant_type',
 'n_affected_exons',
 'variant_filter',
 'gene',
 'decon_correlation',
 'exon_start',
 'exon_end',]]



# Make sure we have at least one numeric column for the group command below
merged_df['start'] = pd.to_numeric(merged_df['start'])

# we need to add a row with NAs if a sample has not had any calls for that caller
# For example if no manta calls add row like sample_name,manta,NA,NA etc.
grouped = merged_df.groupby(['sample_name', 'caller']).mean().reset_index()

decon_calls = {}
manta_calls = {}

# create dictionaries of samples with calls
for row in grouped.itertuples():
    
    if row.caller == 'decon':
        
        decon_calls[row.sample_name] = True
        
    if row.caller == 'manta':
        
        manta_calls[row.sample_name] = True


additional_na = []

# Add a NA row for samples without calls for that caller
for sample in high_cov_sample_list:
    
    if sample not in manta_calls:
        
        additional_na.append([sample, 'manta', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA'])
        
    if sample not in decon_calls:
        
        additional_na.append([sample, 'decon', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA'])
        
additional_df = pd.DataFrame(additional_na, columns=['sample_name',
                                                 'caller',
                                                 'chromosome',
                                                 'start',
                                                 'end',
                                                 'ref',
                                                 'alt',
                                                 'variant_type',
                                                 'n_affected_exons',
                                                 'variant_filter',
                                                 'decon_correlation',
                                                 'exon_start',
                                                 'exon_end'])

# Merge merged df and additional df
final_df = merged_df.append(additional_df, sort=False).sort_values('sample_name')



# Read in gene df
gene_df = pd.read_csv(gene_bed_file, sep='\t', names=['chrom', 'start', 'end', 'gene'])

# Collapse genes df into df with one gene in each row
grouped_gene_df = gene_df.groupby(['chrom', 'gene']).agg({'start': 'min', 'end': 'max'})
grouped_gene_df = grouped_gene_df.reset_index()

# Annotate with gene information
final_df['gene'] = final_df.apply(get_overlapping_genes, axis=1, args=(grouped_gene_df,))

# Add depth information to final df
depth_df = pd.DataFrame()

for depth_file in sample_depths_files:
    
    new_depth_df = pd.read_csv(depth_file, sep='\t')
    depth_df= depth_df.append(new_depth_df)

final_df = final_df.merge(depth_df, how='left', left_on='sample_name', right_on='sample_id')

final_df['sample_depth'] = final_df['mean']

final_df['sample_filter'] = final_df.apply(add_sample_filter, axis=1)

# unify variant type column
final_df['variant_type'] = final_df.apply(fix_variant_type,axis=1)

final_df[['sample_name',
 'caller',
 'chromosome',
 'start',
 'end',
 'ref',
 'variant_type',
 'n_affected_exons',
 'variant_filter',
 'sample_filter',
 'sample_depth',
 'gene',
 'exon_start',
 'exon_end']].to_csv(output, index=False, sep=',')
