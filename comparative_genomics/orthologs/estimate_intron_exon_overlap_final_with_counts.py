from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
import pandas as pd

def read_gff(file_path):
    features = []
    with open(file_path, 'r') as file:
        for line in file:
            if not line.startswith('#'):
                cols = line.strip().split('\t')
                chromosome = cols[0]
                feature_type = cols[2]
                start, end = int(cols[3]), int(cols[4])
                features.append((chromosome, feature_type, FeatureLocation(start, end)))
    return features

def calculate_overlap(chromosome, start, end, features):
    total_exon_length = 0
    total_intron_length = 0
    uncovered_region = 0
    genecount = 0
    exoncount = 0
    introncount = 0

    for feature_chrom, feature_type, location in features:
        if chromosome == feature_chrom and start < location.end and end > location.start:
            overlap_start = max(start, location.start)
            overlap_end = min(end, location.end)
            overlap_length = overlap_end - overlap_start
#            print(f"feature_chrom={feature_chrom}, overlap_start={overlap_start}, feature_type={feature_type}")

            if feature_type == 'exon':
                total_exon_length += overlap_length
                uncovered_region += overlap_length
                exoncount += 1
            elif feature_type == 'intron':
                total_intron_length += overlap_length
                uncovered_region += overlap_length
                introncount += 1
            elif feature_type == 'gene':
                genecount += 1




    # Ensure that uncovered regions are non-negative
#    uncovered_region = max(0, uncovered_region)

    uncovered_region=end - start - uncovered_region
    return total_exon_length, total_intron_length, uncovered_region, genecount, exoncount, introncount

def process_row(row):
    genome1 = row['genome1']
    genome2 = row['genome2']
    start1 = row['startBp1']
    end1 = row['endBp1']
    start2 = row['startBp2']
    end2 = row['endBp2']
    orient = row['orient']
    chromosome1 = row['chr1']
    chromosome2 = row['chr2']

    # Print the relevant information for each row
    #print(f"Processing row: genome1={genome1}, genome2={genome2}, chr1={chromosome1}, chr2={chromosome2}, start1={start1}, end1={end1}, start2={start2}, end2={end2}")

    # Check for missing or invalid values
    if pd.isna(genome1) or pd.isna(genome2) or pd.isna(start1) or pd.isna(end1) or pd.isna(start2) or pd.isna(end2):
        print("Skipping row due to missing or invalid values.")
        return pd.Series({
            'genome1': genome1,
            'genome2': genome2,
            'chr1': chromosome1,
            'chr2': chromosome2,
            'startBp1': start1,
            'endBp1': end1,
            'startBp2': start2,
            'endBp2': end2,
			'orient': orient,
            'exon_length_genome1': 0,
            'intron_length_genome1': 0,
            'exon_length_genome2': 0,
            'intron_length_genome2': 0,
            'uncovered_region_genome1': end1 - start1,
            'uncovered_region_genome2': end2 - start2,
            'genecount_genome1': 0,
            'genecount_genome2': 0,
            'exoncount_genome1': 0,
            'exoncount_genome2': 0,
            'introncount_genome1': 0,
            'introncount_genome2': 0
        })

    gff_file1 = f"{genome1}_ann_w_introns.gff"
    gff_file2 = f"{genome2}_ann_w_introns.gff"

    #print(f"Processing row for genome1={genome1}, genome2={genome2}, chr1={chromosome1}, chr2={chromosome2}, start1={start1}, end1={end1}, start2={start2}, end2={end2}")

    features1 = read_gff(gff_file1)
    features2 = read_gff(gff_file2)

    #printi(f"Reading GFF files: {gff_file1}, {gff_file2}")
    if orient == '-':
        # Swap startBp2 and endBp2
        start2, end2 = end2, start2


    exon_length1, intron_length1, uncovered1, genecount1, exoncount1,introncount1  = calculate_overlap(chromosome1, start1, end1, features1)
    exon_length2, intron_length2, uncovered2, genecount2, exoncount2,introncount2 = calculate_overlap(chromosome2, start2, end2, features2)

#    print(f"Results: exon_length_genome1={exon_length1}, intron_length_genome1={intron_length1}, uncovered_region_genome1={uncovered1}, gene_count1={genecount1}, exoncount1={exoncount1}, introncount1={introncount1}")
#    print(f"Results: exon_length_genome2={exon_length2}, intron_length_genome2={intron_length2}, uncovered_region_genome2={uncovered2}, gene_count2={genecount2}, exoncount2={exoncount2}, introncount2={introncount2}")

    return pd.Series({
        'genome1': genome1,
        'genome2': genome2,
        'chr1': chromosome1,
        'chr2': chromosome2,
        'startBp1': start1,
        'endBp1': end1,
        'startBp2': start2,
        'endBp2': end2,
        'orient': orient,
        'exon_length_genome1': exon_length1,
        'intron_length_genome1': intron_length1,
        'exon_length_genome2': exon_length2,
        'intron_length_genome2': intron_length2,
        'uncovered_region_genome1': uncovered1,
        'uncovered_region_genome2': uncovered2,
        'genecount_genome1': genecount1,
        'genecount_genome2': genecount2,
		'exoncount_genome1': exoncount1,
		'exoncount_genome2': exoncount2,
		'introncount_genome1': introncount1,
		'introncount_genome2': introncount2
    })

# Load your data
#df = pd.read_csv('TEST.csv')
df = pd.read_csv('Cbren_phasedBlks.csv')

# Filter rows where genome1 or genome2 is Cbren
filtered_df = df[(df['genome1'] == 'Cbren') ]

# Apply the processing function to each row
result_df = filtered_df.apply(process_row, axis=1)

#print(result_df)

# Save the result to a new TXT file
#result_df.to_csv('TEST_exon_intron_count.txt', sep='\t', index=False)
result_df.to_csv('GENESPACE_exon_intron_frac_count_Cbren_vs_othersFIX.txt', sep='\t', index=False)
