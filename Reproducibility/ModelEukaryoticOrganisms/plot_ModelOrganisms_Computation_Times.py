
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
import numpy as np
import statistics
import matplotlib.pyplot as plt
import seaborn as sns
import scipy as sp
import pandas as pd
from Bio import SeqIO
from matplotlib.lines import Line2D
import datetime


def main():

    df = {'Computation Time':[],
        'Category':[],
        'Amino Acid':[]}

    proteomes = ['UP000000625_Ecoli_IsoformsIncluded', 'UP000002311_Scerevisiae_IsoformsIncluded', 'UP000001940_Celegans_IsoformsIncluded', 'UP000000803_Dmelanogaster_IsoformsIncluded', 'UP000000437_Drerio_IsoformsIncluded', 'UP000186698_Xlaevis_IsoformsIncluded', 'UP000000589_Mmusculus_IsoformsIncluded', 'UP000005640_Hsapiens_IsoformsIncluded']
    abbrevs = ['E. coli', 'S. cerevisiae', 'C. elegans', 'D. melanogaster', 'D. rerio', 'X. laevis', 'M. musculus', 'H. sapiens']
    index = 0

    h = open('ComputationTimes_ModelOrganisms_IsoformsIncluded.txt')

    for line in h:
        if 'LCD-Composer' in line:
            organism, aa, *remainder = line.rstrip().split('_')
            organism = organism[0] + '. ' + organism[1:]
                
        if 'Start time' in line:
            items = line.rstrip().split(' ')
            start_dt = items[2] + ' ' + items[3]
            start_dtobj = datetime.datetime.strptime(start_dt, '%Y-%m-%d %H:%M:%S.%f')
            
            #PARSE NEXT LINE, WHICH IS END TIME
            items = h.readline().rstrip().split(' ')
            end_dt = items[2] + ' ' + items[3]
            end_dtobj = datetime.datetime.strptime(end_dt, '%Y-%m-%d %H:%M:%S.%f')
            
            delta_t = end_dtobj - start_dtobj
            seconds = delta_t.total_seconds()

            df['Computation Time'].append(seconds)
            df['Amino Acid'].append( aa )
            df['Category'].append( organism )

    index += 1
    h.close()

    new_df = {'Computation Time':[],
        'Category':[],
        'Standard Deviation':[],
        'Amino Acid':[],
        'Proteome Size':[]}
        
    sizes_df, sizes = get_proteome_size()

    df = pd.DataFrame.from_dict(df)

    for i in range(len(proteomes)):
        temp_df = df[df['Category'] == abbrevs[i]]
        times = list(temp_df['Computation Time'])
        if len(times) > 0:
            new_df['Computation Time'].append( sp.nanmean(times) )
        else:
            new_df['Computation Time'].append( np.nan )
        new_df['Standard Deviation'].append( np.nanstd(times) )
        new_df['Category'].append( abbrevs[i] )
        new_df['Amino Acid'].append(aa)
        new_df['Proteome Size'].append(sizes[i])

    plot_proteomes_comptime(new_df)
    plot_scalability(new_df)


def plot_proteomes_comptime(new_df):

    sns.barplot(x=new_df['Category'], y=new_df['Computation Time'], color='#1f77b4', yerr=new_df['Standard Deviation'])
    plt.xticks(fontname='Arial', fontsize=16, rotation=90)
    plt.yticks(fontname='Arial', fontsize=16)
    plt.xlabel('Proteome', fontname='Arial', fontsize=18)
    plt.ylabel('Computation Time\n(seconds)', fontname='Arial', fontsize=18)
    fig = plt.gcf()
    fig.set_size_inches(7, 5)
    plt.savefig('Fig S6B - ModelOrganisms_LCD-Composer_ComputationTimes_IsoformsIncluded.tiff', bbox_inches='tight', dpi=600)
    plt.close()
    
def get_proteome_size():
    proteomes = ['UP000000625_Ecoli_IsoformsIncluded', 'UP000002311_Scerevisiae_IsoformsIncluded', 'UP000001940_Celegans_IsoformsIncluded', 'UP000000803_Dmelanogaster_IsoformsIncluded', 'UP000000437_Drerio_IsoformsIncluded', 'UP000186698_Xlaevis_IsoformsIncluded', 'UP000000589_Mmusculus_IsoformsIncluded', 'UP000005640_Hsapiens_IsoformsIncluded']
    abbrevs = ['E. coli', 'S. cerevisiae', 'C. elegans', 'D. melanogaster', 'D. rerio', 'X. laevis', 'M. musculus', 'H. sapiens']
    df = {}
    sizes = []
    for proteome in proteomes:
        h = open(proteome + '.fasta')
        lengths = []
        num_prots = 0
        for seq_record in SeqIO.parse(h, 'fasta'):
            seq = str(seq_record.seq)
            if seq[-1] == '*':
                seq = seq[:-1]
            lengths.append( len(seq) )
            num_prots += 1

        abbrev = abbrevs[ proteomes.index(proteome) ]
        df[abbrev] = sum(lengths)
        sizes.append(sum(lengths))
        h.close()
        
    return df, sizes
    
def plot_scalability(new_df):
    abbrevs = ['E. coli', 'S. cerevisiae', 'C. elegans', 'D. melanogaster', 'D. rerio', 'X. laevis', 'M. musculus', 'H. sapiens']
    palette = sns.color_palette("colorblind")
    sns.regplot(x=new_df['Proteome Size'], y=new_df['Computation Time'], line_kws={'zorder':0})
    plt.scatter(x=new_df['Proteome Size'], y=new_df['Computation Time'], color=palette[:8], s=80, zorder=1)
    plt.xlabel('Proteome Size\n(# of amino acids, in millions)', fontname='Arial', fontsize=18)
    plt.ylabel('Mean Computation Time\n(seconds)', fontname='Arial', fontsize=18)
    plt.xticks([x for x in range(0, 27000000, 5000000)], labels=[x for x in range(0, 30, 5)], fontname='Arial', fontsize=16)
    plt.yticks(fontname='Arial', fontsize=16)

    legend_elements = []
    for i in range(len(abbrevs)):
        legend_elements.append( Line2D([0], [0], marker='o', color='w', markerfacecolor=palette[i], label='$\it{'+abbrevs[i]+'}$', markersize=10) )

    leg = plt.legend(handles=legend_elements, loc=2, prop={'family':'Arial', 'size':12}, handletextpad=0.0, borderpad=0.4, handlelength=2, labelspacing=0.15)
    leg.set_title('Proteome', prop={'family':'Arial', 'size':14})
    plt.tight_layout()
    fig = plt.gcf()
    fig.set_size_inches(6, 5)
    plt.savefig('Fig S6C - ModelOrganisms_ComputationTimeScalability_IsoformsIncluded.tiff', bbox_inches='tight', dpi=600)
    plt.close()
   
    
if __name__ == '__main__':
    main()
    