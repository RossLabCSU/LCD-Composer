

amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
import numpy as np
from Bio import SeqIO
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch
import pickle
taxonid_to_organism = pickle.load(open('TaxonID_to_Organism.dat', 'rb'))

def main():
    
    life_domain = 'Archaea'

    proteomes = get_proteomes(life_domain)
    total_lcd_contents = {proteome:0 for proteome in proteomes}
    lcd_contents_bytype = {proteome:{aa:0 for aa in amino_acids} for proteome in proteomes}
    
    h = open(life_domain + '_LCDcontent_Statistics.tsv')
    
    header = h.readline()
    header = h.readline()
    for line in h:
        proteome, organism_name, lcd_class, minimum, percentile25, median, percentile75, maximum, average, total_percentage = line.rstrip().split('\t')
        total_percentage = float(total_percentage)
        total_lcd_contents[proteome] += total_percentage
        lcd_contents_bytype[proteome][lcd_class] += total_percentage
    h.close()
        
    plotting(total_lcd_contents, lcd_contents_bytype, proteomes, life_domain)
        
    write_file(total_lcd_contents, lcd_contents_bytype, proteomes)

    
def write_file(total_lcd_contents, lcd_contents_bytype, proteomes):

    output = open('Table S1 - ' + life_domain + ' Percentage_LCDcoverage_wideform.tsv', 'w')
    output.write('\t'.join( ['Proteome', 'Organism'] + list(amino_acids) + ['Cumulative Percentage of Proteome Classified as LCD (NOTE: likely inflated due to overlapping LCDs not being considered in the calculations'] ) + '\n')
    
    for proteome in proteomes:
        #NOTE THAT ISOFORMS IN '_additional' FILES WERE ALREADY INCORPORATED INTO THE PERCENTAGES OUTPUT IN THE '_LCDcontent_Statistics.tsv' FILES. THEREFORE, THEY CAN BE SKIPPED HERE, THOUGH THEY ARE STILL INCLUDED IN THE ACTUAL DATA. THIS IS SIMPLY TO MAKE THE OUTPUT MORE READABLE AND INCLUDE ONLY ONE MENTION OF EACH ORGANISM.
        if '_additional' in proteome:
            continue
        uniprot, taxonid = proteome.split('_')
        organism = taxonid_to_organism[taxonid]
        output.write('\t'.join( [proteome, organism] + [str(lcd_contents_bytype[proteome][aa]) for aa in amino_acids] + [str(total_lcd_contents[proteome])] ) + '\n')
        
    output.close()
    
        
def plotting(total_lcd_contents, lcd_contents_bytype, proteomes, life_domain):

    colors = ['#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c', '#98df8a', '#d62728', '#ff9896', '#9467bd', '#c5b0d5', '#8c564b', '#c49c94', '#e377c2', '#f7b6d2', '#7f7f7f', '#c7c7c7', '#bcbd22', '#dbdb8d', '#17becf', '#9edae65']
    colors = [hex_to_floats(hex) for hex in colors]

    num_proteomes = 30  #NUMBER OF PROTEOMES THAT YOU WISH TO INCLUDE ON THE PLOT

    totals = [total_lcd_contents[proteome] for proteome in proteomes]
    
    totals, sorted_proteomes = zip(*sorted(zip(totals, proteomes), reverse=True))
    
    xlabels = []
    x_index = 0
    for proteome in sorted_proteomes[:num_proteomes]:
        uniprot, taxonid = proteome.split('_')
        common_name = taxonid_to_organism[ taxonid ]
        common_name, *strain_info = common_name.split('(strain')     #HANDLES ORGANISMS WITH LENGTHY NAMES DUE TO SPECIFIC STRAIN INFORMATION IN PARENTHESES
        common_name, *strain_info = common_name.split(' = ')     #HANDLES ORGANISMS WITH LENGTHY NAMES DUE TO SPECIFIC STRAIN INFORMATION IN PARENTHESES
        xlabels.append( common_name )
        
        bottom = 0
        height_sum = 0
        colors_index = 0
        for aa in amino_acids:
            plt.bar(x=x_index, height=lcd_contents_bytype[proteome][aa], bottom=height_sum, color=colors[colors_index])
            colors_index += 1
            height_sum += lcd_contents_bytype[proteome][aa]
            
        x_index += 1

    #COMMENT BACK IN FOR LEGEND=========================================================================================
    # leg_items = [Patch(facecolor=colors[i], label=amino_acids[i]) for i in range(len(amino_acids))]
    
    # # leg = plt.legend([x for x in range(1, 8)], prop={'family':'Arial', 'size':14}, loc=2, bbox_to_anchor=(1.0, 1.0), title='Number of Organisms\nwith GO Enriched')
    # leg = plt.legend(handles=leg_items, prop={'family':'Arial', 'size':10}, loc=2, bbox_to_anchor=(1.0, 1.0), title='LCD Class', handletextpad=0.2)
    # # leg.set_title('Number of Organisms\nwith GO Enriched', fontname='Arial', fontsize=12)
    # leg.set_title('LCD Class')
    # plt.setp(leg.get_title(), multialignment='center', fontname='Arial', fontsize=14)
    #COMMENT BACK IN FOR LEGEND=========================================================================================
    
    plt.xticks([x for x in range(num_proteomes)], labels=xlabels, fontname='Arial', fontsize=10, rotation=90)
    plt.yticks(fontname='Arial', fontsize=14)
    plt.xlabel('Organism', fontname='Arial', fontsize=16)
    plt.ylabel('Percentage of Proteome\nClassified as LCD', fontname='Arial', fontsize=16)
    plt.title(life_domain, fontname='Arial', fontsize=16)
    
    max_total = max(totals)
    plt.ylim(0, max_total+(max_total*0.1))
    plt.xlim(-0.5, num_proteomes-0.5)
    
    plt.tight_layout()
    
    fig = plt.gcf()
    fig.set_size_inches(6, 10)
    plt.savefig('Fig 5A - ' + life_domain + '_Top' + str(num_proteomes) + '_LCDproteins.tiff', bbox_inches='tight', dpi=600)
    plt.close()


def hex_to_floats(hex_code):
    return tuple(int(hex_code[i:i + 2], 16) / 255. for i in (1, 3, 5))
    

def get_proteomes(life_domain):

    h = open(life_domain + '.txt')
    proteomes = []
    for line in h:
        items = line.rstrip().split(' ')
        filename = items[-1]
        if filename.endswith('.fasta.gz') and 'DNA' not in filename:
            accession, junk, extension = filename.split('.')
            proteomes.append( accession )
    h.close()
    
    return proteomes


if __name__ == '__main__':
    main()
