
aa_names = {'A':'Alanine (A)', 'C':'Cysteine (C)', 'D':'Aspartic Acid (D)', 'E':'Glutamic Acid (E)',
        'F':'Phenylalanine (F)', 'G':'Glycine (G)', 'H':'Histidine (H)', 'I':'Isoleucine (I)', 
        'K':'Lysine (K)', 'L':'Leucine (L)', 'M':'Methionine (M)', 'N':'Asparagine (N)',
        'P':'Proline (P)', 'Q':'Glutamine (Q)', 'R':'Arginine (R)', 'S':'Serine (S)',
        'T':'Threonine (T)', 'V':'Valine (V)', 'W':'Tryptophan (W)','Y':'Tyrosine (Y)'}
        
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Patch
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
            
def main():
    
    df = {'Category':[],
        'Amino Acid':[],
        'Composition':[]}
    
    for domain in ['Archaea', 'Bacteria', 'Eukaryota', 'Viruses']:
        df = get_organism_data(df, domain)

    df = pd.DataFrame.from_dict(df)
    stacked_barplot(df)

    
def stacked_barplot(df):

    colors = sns.color_palette('colorblind')
    
    aa_index = 1
    fig, ax = plt.subplots(5,4,sharex=True, sharey=True, figsize=(6.5,10))

    output = open('TableS2_LCDcoverage_StackedBarchart_DATA.csv', 'w')

    for aa in amino_acids:

        ax = plt.subplot(5,4, aa_index)
        small_df = df[df['Amino Acid'] == aa]
        color_index = 0
        x_index = 0
        for domain in ['Viruses', 'Archaea', 'Bacteria', 'Eukaryota']:
            plotting_df = small_df[small_df['Category'] == domain]
            vals = list(plotting_df['Composition'])
            bar_sizes, bar_labels = get_bar_sizes(vals)
            
            if aa == 'A' and domain == 'Viruses':
                output.write(','.join( ['LCD Class', 'Domain of Life'] + bar_labels ) + '\n')
            
            output.write(','.join( [aa, domain] + [str(x) for x in bar_sizes] ) + '\n')
            
            height_sum = 0
            colors_index = 0
            for val in bar_sizes:
                plt.bar(x=x_index, height=val, bottom=height_sum, color=colors[colors_index])
                colors_index += 1
                height_sum += val
            x_index += 1  

        plt.title(aa_names[aa])
        if aa_index in (2,3,4,6,7,8,10,11,12,14,15,16,18,19,20):
            ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
            ax.set_yticklabels([])
        else:
            ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
            ax.set_yticklabels(['0.0', '0.2', '0.4', '0.6', '0.8', '1.0'], fontname='Arial', fontsize=12)
        if aa_index in range(1,17):
            ax.set_xticks([0,1,2,3])
            ax.set_xticklabels([])
        if aa_index in range(17,21):
            ax.set_xticks([0,1,2,3])
            ax.set_xticklabels(['Viruses', 'Archaea', 'Bacteria', 'Eukaryota'], rotation=90, fontname='Arial', fontsize=12)
        
        aa_index += 1

    plt.subplots_adjust(wspace=-0.1, hspace=0)
        
    #LEGEND
    legend_elements = [Patch(facecolor=colors[i], label=bar_labels[i]) for i in range(len(bar_labels))]
    leg = fig.legend(handles=legend_elements, loc="center left", bbox_to_anchor=(1,0.5), handletextpad=0.2, title='LCD Coverage\n(Percentage of Proteome)', prop={'family':'Arial', 'size':10})
    leg.set_title('LCD Content (%)',prop={'size':12, 'family':'Arial'})

    fig.text(-0.01, 0.5, 'Proportion of Proteomes in LCD Content Bin', va='center', rotation='vertical', fontname='Arial', fontsize=14)
    plt.tight_layout()
    plt.savefig('Fig 4 - LCDcoverage_StackedBarcharts.tiff', bbox_inches = 'tight', dpi=600)
    plt.close() 
        
        
def get_bar_sizes(vals):
    
    #INITIALIZE THE BAR SIZES LIST BY COUNTING THE NUMBER OF PROTEOMES FOR WHICH 0% OF THE PROTEOME WAS CLASSIFIED AS LCD
    bar_sizes = [vals.count(float('0')) / len(vals)]
    bar_labels = ['None (0%)', 'Extremely Low (0-0.1%)', 'Very Low (0.1-0.5%)', 'Low (0.5-2%)', 'Medium (2-5%)', 'High (5-10%)', 'Very High (10-15%)', 'Extremely High (>15%)']
    for bounds in [(0, 0.1), (0.1,0.5), (0.5,2), (2,5), (5,10), (10,15), (15,100)]:
        size = sum([1 if val > bounds[0] and val <= bounds[1] else 0 for val in vals]) / len(vals)
        bar_sizes.append(size)
        
    return bar_sizes, bar_labels
    
    
def get_organism_data(df, domain):
    
    h = open(domain + '_LCDcontent_Statistics.tsv')
    
    header = h.readline()
    header = h.readline()
    
    for line in h:
        items = line.rstrip().split('\t')
        lcd_class = items[2]
        composition = float(items[-1])
        
        df['Category'].append(domain)
        df['Amino Acid'].append(lcd_class)
        df['Composition'].append(composition)
        
    h.close()
    
    return df
     

            
if __name__ == '__main__':
    main()
    