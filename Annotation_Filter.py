## SETUP
# Packages
import os
import pandas as pd
import subprocess

# User to change - paths
out_prefix = "seq_1"
root_path = "/home/"
amr_finder_path = root_path + "ARG_Annotation/"+ out_prefix + "_ARG.txt" # path to AMRFinderPlus tabular output file
fasta_path = root_path + out_prefix + ".fasta" # .fasta for sequence to be annotated
mge_annotation_outpath = root_path + "MGE_Annotation/" # output directory
mgedb_fasta_path = root_path + "MGE_db_dd.fasta" # path to mge database (.fasta)
additional_elements_path = root_path + "Additional_Elements/" # path to directory with .fasta files for elements not in mge database

# Directories
os.makedirs(mge_annotation_outpath, exist_ok = True)
os.chdir(mge_annotation_outpath)

print("Aligning input sequence to MGEs database")

# Getting ARG Information
amr_finder = pd.read_table(amr_finder_path, usecols=[2,3,4,5,6,10,11,13,14,15,16,17])

# Running MGE BLAST and importing results
# Combining MGEdb with the additional elements folder and running BLAST
subprocess.run( f"cat {mgedb_fasta_path} {additional_elements_path}* >> mge_db_temp.fasta", shell=True)
subprocess.run( f"makeblastdb -dbtype nucl -in mge_db_temp.fasta -parse_seqids -out mge_db_temp_blastdb; rm mge_db_temp.fasta", shell=True)
subprocess.run( f"blastn -task blastn -query {fasta_path} -db mge_db_temp_blastdb -perc_identity 90 -max_target_seqs 7000 -outfmt '6 std qlen slen' -out {out_prefix}_MGE.blastout", shell=True)
subprocess.run( f"rm mge_db_temp_blastdb*", shell=True)
blast_field_names = ['MGE', 'Identity', 'Alignment_Length', 'Acc_Start', 'Acc_End', 'MGE_Start', 
                    'MGE_End', 'Acc_Length', 'MGE_Length']
mge_blast_out = pd.read_table(mge_annotation_outpath + out_prefix + "_MGE.blastout", header=None, 
                            usecols=[1,2,3,6,7,8,9,12,13], names = blast_field_names)
# Adding a coverage column to the MGE output
mge_blast_out['MGE_Coverage'] = mge_blast_out['Alignment_Length'] / mge_blast_out['MGE_Length']
print("Originally " + str(len(mge_blast_out)) + " alignments")


print("\nFiltering alignments to MGE database")

## STAGE 1: IS26 CHECK

# Setup
mge_blast_out_filt1 = mge_blast_out
mge_blast_out_keep = pd.DataFrame()

# Take all IS26-like elements > 99% identity
IS26_like_df_raw = pd.concat([mge_blast_out[mge_blast_out['MGE'].str.contains("IS26|IS15")], 
                        mge_blast_out[mge_blast_out['MGE'].str.contains(r'^(?=.*IS26)(?=.*IS15)')]]
                        ).drop_duplicates().sort_values(by = ['Acc_Start'], ignore_index = True)
IS26_like_df = IS26_like_df_raw.loc[IS26_like_df_raw['Identity'] > 99]

# If there are IS26s present
coordinates = [0,0]

if (len(IS26_like_df)) > 0:
    for index, complete_IS26 in IS26_like_df.iterrows():
        # If the current IS26 doesn't have coordinates within 20 bp of the previous one
        if not (complete_IS26['Acc_Start'] > min(coordinates) - 20) & \
            (complete_IS26['Acc_Start'] < max(coordinates) + 20) & \
            (complete_IS26['Acc_End'] > min(coordinates) - 20) & \
            (complete_IS26['Acc_End'] < max(coordinates) + 20):

            # Identify all annotations with coordinates within 20 bp of the current IS26 
            coordinates = [complete_IS26['Acc_Start'], complete_IS26['Acc_End']]
            within_current = mge_blast_out.loc[((mge_blast_out['Acc_Start'].between(min(coordinates) - 20, max(coordinates) + 20)) &
                                                (mge_blast_out['Acc_End'].between(min(coordinates) - 20, max(coordinates) + 20)))]
            
            # Pick the best IS26 in that region (first by identity, then coverage) and add to keep pile
            current_region_IS26s = pd.concat([within_current[within_current['MGE'].str.contains("IS26|IS15")], 
                                        within_current[within_current['MGE'].str.contains(r'^(?=.*IS26)(?=.*IS15)')]]
                                        ).drop_duplicates().sort_values(by = 'Identity', ascending = False)

            max_IS26_id = max(current_region_IS26s.loc[:]['Identity'].tolist())
            best_IS26 = current_region_IS26s.loc[(current_region_IS26s['Identity'] == max_IS26_id)]            
            if len(best_IS26) > 1:
                max_cov = max(best_IS26.loc[:]['MGE_Coverage'].tolist())
                best_IS26 = best_IS26.loc[best_IS26['MGE_Coverage'] == max_cov]

            mge_blast_out_keep = pd.concat([mge_blast_out_keep, best_IS26])

            # Filter out alignments to the same region of query sequence as current IS26 (+/- 20 bp)
            IS26_coordinates = [best_IS26.iloc[0]['Acc_Start'], best_IS26.iloc[0]['Acc_End']]
            mge_blast_out_filt1 = mge_blast_out_filt1.loc[((mge_blast_out_filt1['Acc_Start'] < min(IS26_coordinates) - 20) | 
                                    (mge_blast_out_filt1['Acc_End'] > max(IS26_coordinates) + 20)) | 
                                    ((mge_blast_out_filt1['Acc_End'] < min(IS26_coordinates) - 20) | 
                                    (mge_blast_out_filt1['Acc_Start'] > max(IS26_coordinates) + 20))]
else:
    # If there are no IS26s, move along
    mge_blast_out_filt1 = mge_blast_out

print("After filtering for same coordinates as IS26s, " + str(len(mge_blast_out_filt1)) + " remaining")



## STAGE 2: FULL LENGTH & 100% IDENTITY KEPT

# Everything full length and 100% in the keep pile
mge_blast_out_full = mge_blast_out_filt1.loc[(mge_blast_out_filt1['Identity'] == 100) & (mge_blast_out_filt1['Alignment_Length'] == 
                mge_blast_out_filt1['MGE_Length'])]
mge_blast_out_keep2 = pd.concat([mge_blast_out_keep, mge_blast_out_full])

# Filter out alignments to the same region of query sequence as keep pile items (+/-20 bp)
mge_blast_out_filt2 = mge_blast_out_filt1
for index, kept_element in mge_blast_out_keep2.iterrows():
    coordinates = [kept_element['Acc_Start'], kept_element['Acc_End']]
    mge_blast_out_filt2 = mge_blast_out_filt2.loc[((
        mge_blast_out_filt2['Acc_Start'] < min(coordinates) - 20) | 
        (mge_blast_out_filt2['Acc_End'] > min(coordinates) + 20)) | 
        ((mge_blast_out_filt2['Acc_End'] < max(coordinates) - 20) | 
        (mge_blast_out_filt2['Acc_Start'] > max(coordinates) + 20))]

print("After filtering for same coordinates as perfect alignments, " + str(len(mge_blast_out_filt2)) + " remaining")




## STAGE 3: REMOVAL OF ANY ALIGNMENTS JUST TO RESISTANCE GENES

mge_blast_out_filt3 = mge_blast_out_filt2

# Remove alignments to the same region of query sequence as ARGs (+/-20 bp)
for index, ARG in amr_finder.iterrows():
    coordinates = [ARG['Start'], ARG['Stop']]
    mge_blast_out_filt3 = mge_blast_out_filt3.loc[((mge_blast_out_filt3['Acc_Start'] < min(coordinates) - 20) | 
                                (mge_blast_out_filt3['Acc_End'] > min(coordinates) + 20)) | 
                                ((mge_blast_out_filt3['Acc_End'] < max(coordinates) - 20) | 
                                (mge_blast_out_filt3['Acc_Start'] > max(coordinates) + 20))]

print("After filtering for same coordinates as ARGs, " + str(len(mge_blast_out_filt3)) + " remaining")




## STAGE 4: REMOVAL OF ANYTHING WITH < 10% COVERAGE

mge_blast_out_filt4 = mge_blast_out_filt3.loc[mge_blast_out_filt3['MGE_Coverage'] >= 0.1]

print("After filtering out elements with <10% coverage, " + str(len(mge_blast_out_filt4)) + " remaining")




## STAGE 5: TAKE LONGEST ALIGNMENT AT EACH LOCATION, PLUS ANY INTERNALS WITH >= LONGEST'S IDENTITY 
## OR >=60% COVERAGE

# Iterative
# Starting with the longest alignment, take every alignment that has query sequence coordinates in 
# within the longest (+/- 20bp), filter out anything lower identity than longest or <60% coverage
mge_blast_out_filt4_sorted = pd.concat([mge_blast_out_filt4, mge_blast_out_keep2]).sort_values(
                                by = ['Alignment_Length'], ascending = False, ignore_index = True)
mge_blast_out_filt5 = pd.DataFrame()

# Until everything has either been kept in _filt5 or removed, 
while len(mge_blast_out_filt4_sorted) > 0:

    # Get coordinates and identity of longest
    coordinates = [mge_blast_out_filt4_sorted.head(n = 1)['Acc_Start'].tolist()[0], 
                    mge_blast_out_filt4_sorted.head(n = 1).iloc[:]['Acc_End'].tolist()[0]]
    current_id = mge_blast_out_filt4_sorted.head(n = 1).iloc[:]['Identity'].tolist()[0]

    # Keep anything internal to longest (+/- 20 bp) with higher identity or >=60% coverage
    kept_internals = mge_blast_out_filt4_sorted.loc[
        (mge_blast_out_filt4_sorted['Acc_Start'].between(min(coordinates) - 20, max(coordinates) + 20)) &
        (mge_blast_out_filt4_sorted['Acc_End'].between(min(coordinates) - 20, max(coordinates) + 20)) & 
        ((mge_blast_out_filt4_sorted['Identity'] >= current_id) | (mge_blast_out_filt4_sorted['MGE_Coverage'] >= 0.6))]
    

    # Regenerate remaining alignment dataframe without alignments in the region that was just filtered
    mge_blast_out_filt4_sorted = mge_blast_out_filt4_sorted.loc[(
        (mge_blast_out_filt4_sorted['Acc_Start'] <= min(coordinates) - 20) | 
        (mge_blast_out_filt4_sorted['Acc_Start'] >= max(coordinates) + 20)) | (
        (mge_blast_out_filt4_sorted['Acc_End'] <= min(coordinates) - 20) | 
        (mge_blast_out_filt4_sorted['Acc_End'] >= max(coordinates) + 20))]
    
    # Add retained alignments (longest + it's kept internals) to new filtered dataframe (filt5)
    mge_blast_out_filt5 = pd.concat([mge_blast_out_filt5, kept_internals])

print("After filtering for internals of longest elements, " + str(len(mge_blast_out_filt5)) + " remaining")





## STAGE 6: COMBINE THE KEEP AND MAYBE PILES FOR A FINAL LIST OF PLAUSIBLE ANNOTATIONS
final_filtered_elements = pd.concat([mge_blast_out_filt5, mge_blast_out_keep2]).drop_duplicates()
print("After all automated filtering, " + str(len(mge_blast_out_filt5)) + " remaining")




## SAVE: FILTERED DATAFRAME TO .TSV
outpath = mge_annotation_outpath + out_prefix + "_MGE_AutoFilt.tsv"
final_filtered_elements.to_csv(outpath, sep='\t', header=False, index=False)
