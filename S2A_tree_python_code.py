import pandas as pd
import numpy as np
from Bio import Phylo
from collections import Counter


gtdb_folder = '/scr/user/eryixian/tree_skin/'
tree_file = gtdb_folder+'S2A_tree_coded.nwk'
tree = Phylo.read(tree_file, "newick")



excel_df = pd.read_csv(gtdb_folder+'Fig_S2A_metadata.csv', index_col=0,sep=',')
excel_df['phylum'] = [i.split(';p__')[1].split(';')[0] if ';p' in i else i for i in excel_df['classification']]##
novelty_list = excel_df['novel']


from collections import Counter, OrderedDict
import seaborn as sns

phylum_counts = Counter(excel_df['phylum'])
ordered_phylum = [it[0] for it in OrderedDict(sorted(phylum_counts.items(),key=lambda t:t[1])[::-1]).items() 
                  if it[0] != ""]
colors = sns.color_palette("Paired").as_hex()+sns.color_palette("Set2").as_hex()
phylum_color = {p:colors[i] for i, p in enumerate(ordered_phylum)}
phylum_color

phyla_list = []
with open(tree_file.replace('.nwk', '.annotation.color_strip.txt'), 'w') as wf:
    wf.write('DATASET_COLORSTRIP\n')
    wf.write('SEPARATOR SPACE\n')
    wf.write('DATASET_LABEL color_strip\n')
    wf.write('COLOR #ff0000\n')
    wf.write('COLOR_BRANCHES 0\n')
    wf.write('LEGEND_TITLE Phylum\n')
    wf.write('LEGEND_POSITION_X 100\n')
    wf.write('LEGEND_POSITION_Y 100\n')
    wf.write('LEGEND_SHAPES ' + ' '.join(['1'] * len(phylum_color)) + '\n')
    wf.write('LEGEND_COLORS')
    for p in phylum_color:
        wf.write(' ' + phylum_color[p])
    wf.write('\nLEGEND_LABELS')
    for p in phylum_color:
        wf.write(' ' + p)
    wf.write('\nLEGEND_SHAPE_SCALES ' + ' '.join(['1'] * len(phylum_color)) + '\n')
    wf.write('STRIP_WIDTH 25\n')
    wf.write('MARGIN 0\n')
    wf.write('BORDER_WIDTH 0\n')
    wf.write('MARGIN 0\n')
    wf.write('DATA\n')
    for name in excel_df.index:
        phylum = excel_df.loc[name, 'phylum']
        color = phylum_color.get(phylum, '#000000')  
        wf.write(f'{name} {color} {phylum}\n')
        phyla_list.append(phylum)



phylum_counts = Counter(excel_df['novel'])
ordered_phylum = [it[0] for it in OrderedDict(sorted(phylum_counts.items(),key=lambda t:t[1])[::-1]).items() 
                  if it[0] != ""]
colors = ["#ff0000", "#fdefb2"]
phylum_color = {p:colors[i] for i, p in enumerate(ordered_phylum)}
phylum_color


phyla_list = []
with open(tree_file.replace('.nwk', '.annotation.color_strip2.txt'), 'w') as wf:
    wf.write('DATASET_COLORSTRIP\n')
    wf.write('SEPARATOR SPACE\n')
    wf.write('DATASET_LABEL color_strip\n')
    wf.write('COLOR #ff0000\n')
    wf.write('COLOR_BRANCHES 0\n')
    wf.write('LEGEND_TITLE Source\n')
    wf.write('LEGEND_POSITION_X 100\n')
    wf.write('LEGEND_POSITION_Y 100\n')
    wf.write('LEGEND_SHAPES ' + ' '.join(['1'] * len(phylum_color)) + '\n')
    wf.write('LEGEND_COLORS')
    for p in phylum_color:
        wf.write(' ' + phylum_color[p])
    wf.write('\nLEGEND_LABELS')
    for p in phylum_color:
        wf.write(' ' + p)
    wf.write('\nLEGEND_SHAPE_SCALES ' + ' '.join(['1'] * len(phylum_color)) + '\n')
    wf.write('STRIP_WIDTH 25\n')
    wf.write('MARGIN 0\n')
    wf.write('BORDER_WIDTH 0\n')
    wf.write('MARGIN 0\n')
    wf.write('DATA\n')
    for name in excel_df.index:
        novel = excel_df.loc[name, 'novel']  
        color = phylum_color.get(novel, '#000000') 
        wf.write(f'{name} {color} {novel}\n')
        phyla_list.append(novel)




output_file = tree_file.replace('.nwk', '.annotation.simplebar.txt')

with open(output_file, 'w') as wf:
    wf.write('DATASET_SIMPLEBAR\n')
    wf.write('SEPARATOR COMMA\n')
    wf.write('DATASET_LABEL,genome counts\n')
    wf.write('COLOR,#808080\n')
    wf.write('WIDTH,200\n')
    wf.write('MARGIN,0\n')
    wf.write('HEIGHT_FACTOR,1\n')
    wf.write('BAR_SHIFT,0\n')
    wf.write('BAR_ZERO,0\n')
    wf.write('DATA\n')

    # Write bar data
    for genome_id, row in excel_df.iterrows():
        cleaned_id = genome_id.replace('+', '_')  
        wf.write(f"{cleaned_id},{row['count']}\n")
