import pandas as pd
import sys
import os



def load_orientation_data(gene_name, file_path):
    data = pd.read_csv(file_path, sep='\t', header=None,
                       names=['intron_size', 'orientation', 'intron_index', 'gene_name'])
#    print(data[data['gene_name'] == gene_name])
    return data[data['gene_name'].isin([gene_name])]

def load_phase_data(gene_name, file_path):
    data = pd.read_csv(file_path, sep='\t', header=None,
                       names=['phase', 'cds_name', 'gene_name'])
    return data[data['gene_name'].isin([gene_name])]

def load_splicing_data(gene_name, file_path):
    data = pd.read_csv(file_path, sep=' ', header=None,
                       names=['splicing_site', 'gene_name', 'donor', 'acceptor'])
    #print(data.iloc[0])
    return data[data['gene_name'].isin([gene_name])]

def load_gene_data(file_path):
    data = pd.read_csv(file_path, sep='\t', header=None,
			names=['species_name', 'gene_name', 'gene_structure'],dtype={'gene_structure': str})

    # Extract the first part of the input file as a variable
    first_part_variable = os.path.splitext(os.path.basename(file_path))[0]

    return data, first_part_variable

def find_orientation_entry(gene_name, orientation_data):
    return orientation_data[orientation_data['gene_name'].isin([gene_name])].iloc[0]

def get_intron_sizes(gene_name, orientation_entry):
    return orientation_entry['intron_size'].tolist()

def get_intron_sizes_reverse(gene_name, orientation_entry):
    return orientation_entry['intron_size'].tolist()[::-1]

def get_phases(gene_name, phase_data):
    return phase_data[phase_data['gene_name'].isin([gene_name])]['phase'].tolist()

def get_phases_reverse(gene_name, phase_data):
    return phase_data[phase_data['gene_name'].isin([gene_name])]['phase'].tolist()[::-1]

def get_splicing_sites(gene_name, splicing_data, orientation):
    splicing_entries = splicing_data[splicing_data['gene_name'].isin([gene_name])]
    sites = [f"{entry['donor']}-{entry['acceptor']}" for _, entry in splicing_entries.iterrows()]

    return sites[::-1] if orientation == '-' else sites




def replace_ones_with_values(the_list, replacement_values):
    new_list = the_list.copy()

    if '1' in the_list and replacement_values:
        for i, val in enumerate(new_list):
            if val == '1':
                new_list[i] = str(replacement_values.pop(0))
            else:
                new_list[i] = '0'
    return new_list





gene_file = sys.argv[1]

# Load gene data to extract species_name and first part variable
gene_data, first_part_variable = load_gene_data(gene_file)

gene_result_df = pd.DataFrame()

for _, gene_entry in gene_data.iterrows():
        gene_name = gene_entry['gene_name']
        species_name = gene_entry['species_name']
        gene_structure_list = list(str(gene_entry['gene_structure']))
        #print(gene_structure_list)

        orientation_file = f'/projects/phillipslab/ateterina/Cbren/phylogenomics/worm23/genespace_newref/orthofinder/Results_Oct28/Orthogroups/{species_name}_1-to-1_intlen.tab'
        phase_file = f'/projects/phillipslab/ateterina/Cbren/phylogenomics/worm23/genespace_newref/orthofinder/Results_Oct28/Orthogroups/{species_name}_1-to-1_PHASE.tab'
        splicing_file = f'/projects/phillipslab/ateterina/Cbren/phylogenomics/worm23/genespace_newref/orthofinder/Results_Oct28/Orthogroups/{species_name}_1-to-1_5-3combo.ss.tab'

        print("Gene Name:", gene_name)
        orientation_data = load_orientation_data(gene_name, orientation_file)
        phase_data = load_phase_data(gene_name, phase_file)
        splicing_data = load_splicing_data(gene_name, splicing_file)
        orientation_entry = find_orientation_entry(gene_name, orientation_data)
        orientation = orientation_entry['orientation']

        print(orientation)

        if orientation == '+':
            intron_sizes = get_intron_sizes(gene_name, orientation_data)
            phases = get_phases(gene_name, phase_data)
            splicing_sites = get_splicing_sites(gene_name, splicing_data, orientation)
        elif orientation == '-':
            intron_sizes = get_intron_sizes_reverse(gene_name, orientation_data)
            phases = get_phases_reverse(gene_name, phase_data)
            splicing_sites = get_splicing_sites(gene_name, splicing_data, orientation)
        intron_val=replace_ones_with_values(gene_structure_list,intron_sizes)
        phases_val=replace_ones_with_values(gene_structure_list,phases)
        splicing_val=replace_ones_with_values(gene_structure_list,splicing_sites)
        species_df = pd.DataFrame({
			    "gene_structure": gene_structure_list,
				    "intron_sizes": intron_val,
					    "phases": phases_val,
						    "splicing_sites": splicing_val
        })
        gene_result_df = pd.concat([gene_result_df, species_df], axis=1)






# Add the first part variable column to the result DataFrame
gene_result_df['first_part_variable'] = first_part_variable


#print(gene_result_df)
output_file = '7sp_all_intron_info_ALL.FIN.tab'
gene_result_df.to_csv(output_file, index=False, header=False, mode='a',sep="\t")
