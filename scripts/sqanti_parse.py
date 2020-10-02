import pandas as pd
import glob

columns_name = ['sample', 'novel_genes', 'annotated_genes', 
           'full-splice_match', 'novel_not_in_catalog',
           'intergenic', 'novel_in_catalog', 'antisense',
           'genic', 'fusion', 'incomplete-splice_match', 'genic_intron',
           'known_canonical', 'novel_canonical', 'known_non_canonical', 'novel_non_canonical']

output = pd.DataFrame(columns=columns_name)

for filename in list(snakemake.params.path):
    # Classification
    classification_tsv = pd.read_csv(f"{filename}_classification.txt", sep='\t',
                                    dtype={"chrom": "str"})
    
    line = {}
    line['sample'] = filename.split("/")[-1].split('_')[0]
    line['novel_genes'] = classification_tsv['associated_transcript'].value_counts()['novel']
    line['annotated_genes'] = len(classification_tsv['associated_transcript']) - line['novel_genes']
    
    trans_char = classification_tsv['structural_category'].value_counts().to_dict()
    for char in trans_char:
        line[char] = trans_char[char]

    # Junctions
    j_tsv = pd.read_csv(f"{filename}_junctions.txt", sep="\t", dtype={"chrom": "str"})
    j_tsv["junctions_combined"] = j_tsv['junction_category'] + "_" + j_tsv['canonical']
    junctions = j_tsv["junctions_combined"].value_counts().to_dict()
    for junction in junctions:
        line[junction] = junctions[junction]
    
    output = output.append(line, ignore_index=True)

output = output.fillna(0)
output.to_csv(str(snakemake.output.summary), sep='\t', index=False)
print("SQANTI Summary TSV successfully created!")

