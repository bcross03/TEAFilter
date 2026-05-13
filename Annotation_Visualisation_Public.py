## Takes output from Annotation_Filter_Public.py and generates an image representation of filtered alignments

def annotation_visualisation(seq_label, mge_annotation_outpath, mge_path, arg_path, fasta_path, 
                             integron_path, more_height):

    # Loading data
    import pandas as pd
    import os
    from dna_features_viewer import GraphicFeature, GraphicRecord
    import matplotlib.pyplot as plt
    import re
    from Bio import SeqIO

    os.chdir(mge_annotation_outpath)

    # Annotation_Filter_Public.py outputs a reduced number of fields, columns do not match standard BLAST
    mge_colnames = ['el_name', 'identity', 'aln_length', 'acc_st', 'acc_end', 'mge_st', 'mge_end', 
                    'acc_len', 'mge_len']
    col_nums = [0,1,2,3,4,5,6,7,8]

    # Read in MGE and ARG outputs
    mge_blast_out = pd.read_table(mge_path, usecols = col_nums, names = mge_colnames, header=None)
    arg_amr_finder = pd.read_table(arg_path)

    # Setting up plot parameters (heights etc.)
    for record in SeqIO.parse(fasta_path, "fasta"):
        sequence_length = len(record.seq)
    if len(integron_path) > 0:
        integron_output = open(integron_path,"r")
        if re.search("No Integron found", integron_output.read()):
            integron_height = 0.2
        else:
            integron_height = 0.8
    else:
        integron_height = 0.2
    if more_height != False:
        height_div = more_height
    else:
        height_div = 1.8
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(sequence_length / 1000, len(mge_blast_out) + 1 / height_div), 
                                gridspec_kw={"height_ratios": [4, integron_height, 1]})
    fig.tight_layout()


    ### MGEs
    # Creating start, end, strand and label lists
    # Strand
    mge_strands = {}

    for row in range(len(mge_blast_out)):
        if mge_blast_out.loc[row]['mge_st'] < mge_blast_out.loc[row]['mge_end']:
            mge_strands[row] = +1
        elif mge_blast_out.loc[row]['mge_st'] > mge_blast_out.loc[row]['mge_end']:
            mge_strands[row] = -1

    #Start
    mge_starts = mge_blast_out.loc[:]['acc_st'].tolist()

    #End
    mge_ends = mge_blast_out.loc[:]['acc_end'].tolist()

    #Labels
    counter = 0
    mge_labels = []
    for element_name in mge_blast_out.loc[:]['el_name'].tolist():
        identity = str(round(mge_blast_out.loc[counter]['identity'], 1)) + "%"
        coverage = str(round(mge_blast_out.loc[counter]['aln_length'] / mge_blast_out.loc[counter][
                                                                    'mge_len'] * 100, 1)) + "%"
        if min([mge_blast_out.loc[counter]['mge_st'], mge_blast_out.loc[counter]['mge_end']]) == 1:
            left = "L"
        else: 
            left = ""
        if max([mge_blast_out.loc[counter]['mge_st'], mge_blast_out.loc[counter]['mge_end']]) \
                        == mge_blast_out.loc[counter]['mge_len']:
            right = ", R"
        else:
            right = ""
        mge_labels.append(str(element_name) + ", Id " + str(identity) + ", Cov " + str(coverage) + ", Ends " 
                        + str(left) + str(right))
        counter = counter + 1

    # Calculating % Coverage to determine colour
    mge_coverages = []
    for i in range(len(mge_blast_out)):
        mge_coverages.append(round(abs(mge_blast_out.loc[i]['mge_end'] - mge_blast_out.loc[i]['mge_st']) / 
                            mge_blast_out.loc[i]['mge_len'], 3))

    mge_colours = []
    for i in range(len(mge_coverages)):
        if mge_coverages[i] > 0.8:
            mge_colours.append("#A93226")
        elif mge_coverages[i] > 0.6:
            mge_colours.append("#E34234")
        elif mge_coverages[i] > 0.4:
            mge_colours.append("#E67E22")
        elif mge_coverages[i] > 0.2:
            mge_colours.append("#F1C40F")
        else:
            mge_colours.append("#F9E79F")

    # Creating MGE feature list
    mge_features = []

    for i in range(len(mge_starts)):
        mge_features.append(GraphicFeature(start=int(mge_starts[i]), end=int(mge_ends[i]), strand=mge_strands[i], 
                                    color=mge_colours[i], label=mge_labels[i]))
        



    ### ARGs

    # Creating start, end, strand and label lists
    # Strand
    arg_strands = []
    for strand in arg_amr_finder.loc[:]['Strand'].tolist():
        if strand == '+':
            arg_strands.append(1)
        elif strand == '-':
            arg_strands.append(-1)

    #Start
    arg_starts = arg_amr_finder.loc[:]['Start'].tolist()

    #End
    arg_ends = arg_amr_finder.loc[:]['Stop'].tolist()

    #Labels
    counter = 0
    arg_labels = []
    for gene_name in arg_amr_finder.loc[:]['Element symbol'].tolist(): 
    # for gene_name in arg_amr_finder.loc[:]['Gene symbol'].tolist(): # Some AMRFinderPlus versions use this label
        coverage = " (Cov: " + str(arg_amr_finder.loc[counter]['% Coverage of reference']) + "%)"
        # coverage = " (Cov: " + str(arg_amr_finder.loc[counter]['% Coverage of reference sequence']) + "%)" # Some AMRFinderPlus versions use this label
        arg_labels.append(str(gene_name) + str(coverage))
        counter = counter + 1

    # Calculating % Coverage to determine colour
    arg_coverages = []
    for i in range(len(arg_amr_finder)):
        arg_coverages.append(float(arg_amr_finder.loc[i]['% Coverage of reference']))
        # arg_coverages.append(float(arg_amr_finder.loc[i]['% Coverage of reference sequence'])) # Some AMRFinderPlus versions use this label

    arg_colours = []
    for i in range(len(arg_coverages)):
        if arg_coverages[i] > 80:
            arg_colours.append("#A93226")
        elif arg_coverages[i] > 60:
            arg_colours.append("#E74C3C")
        elif arg_coverages[i] > 40:
            arg_colours.append("#E67E22")
        elif arg_coverages[i] > 20:
            arg_colours.append("#F1C40F")
        else:
            arg_colours.append("#F9E79F")

    # Creating ARG feature list
    arg_features = []

    for i in range(len(arg_starts)):
        arg_features.append(GraphicFeature(start=int(arg_starts[i]), end=int(arg_ends[i]), strand=arg_strands[i], 
                                    color=arg_colours[i], label=arg_labels[i]))




    ### INTEGRONS FEATURES (optional)
    if len(integron_path) > 0:
        integron_output = open(integron_path,"r")
        if re.search("No Integron found", integron_output.read()):
            int_features = []
        else:
            integron_finder = pd.read_table(integron_path, comment = "#")

            # Creating start, end, strand and label lists
            # Strand
            int_strands = integron_finder['strand'].tolist()

            #Start
            int_starts = integron_finder['pos_beg'].tolist()

            #End
            int_ends = integron_finder['pos_end'].tolist()

            #Labels and colours
            counter = 0
            int_labels = []
            int_colours = []

            for gene_name in integron_finder.loc[:]['annotation'].tolist():
                feature_type = integron_finder.loc[counter]['type']
                
                #Colours
                if feature_type == "CALIN":
                    int_colours.append("#6CB4EE")
                elif feature_type == "complete":
                    int_colours.append("#00308F")
                elif feature_type == "In0":
                    int_colours.append("#5072A7")

                #Labels
                int_labels.append(str(gene_name))
                counter = counter + 1

            # Creating feature list
            int_features = []

            for i in range(len(int_starts)):
                int_features.append(GraphicFeature(start=int(int_starts[i]), end=int(int_ends[i]), strand=int_strands[i], 
                                            color=int_colours[i], label=int_labels[i]))
    else:
        int_features = []




    # Plotting MGEs
    mge_record = GraphicRecord(sequence_length = sequence_length, features=mge_features)
    mge_record.plot(ax=ax1, figure_width = sequence_length)

    # Plotting Integrons
    mge_record = GraphicRecord(sequence_length = sequence_length, features=int_features)
    mge_record.plot(ax=ax2, figure_width = sequence_length)

    # Plotting ARGS
    arg_record = GraphicRecord(sequence_length = sequence_length, features=arg_features)
    arg_record.plot(ax=ax3, figure_width = sequence_length)


    # Saving image
    save_path = mge_annotation_outpath + seq_label + "_MGE_Viz.png"
    fig.figure.savefig(save_path, bbox_inches='tight')
    
    return(mge_features, arg_features, int_features)





### Setting paths etc. and running visualisation function
root_path = "/home/"
# Name for sequence - will be outfile prefix
seq_label = "seq_1"
# Output location for image
mge_annotation_outpath = root_path + "MGE_Annotation/"
# Path to fasta file for sequence that was annotated
fasta_path = root_path + seq_label + ".fasta"
# Path to output file from Annotation_Filter_Public.py (ends with AutoFilt.tsv)
mge_path = mge_annotation_outpath + seq_label + "_MGE_AutoFilt.tsv"
# Path to AMRFinderPlus output
arg_path = root_path + "ARG_Annotation/"+ seq_label + "_ARG.txt"
# Path to IntegronFinder output (Optional, if not using, value is "")
integron_path = root_path + "Integron_Annotations/Results_Integron_Finder_" \
                    + seq_label + "/" + seq_label + ".integrons"
# Default figure height is pre-set. 
    # If the output is squished, change this to a float smaller than 1.8. 
more_height = False

## Create Image
mge_features, arg_features, int_features = annotation_visualisation(seq_label, 
                                            mge_annotation_outpath, mge_path, arg_path, 
                                            fasta_path, integron_path, more_height)