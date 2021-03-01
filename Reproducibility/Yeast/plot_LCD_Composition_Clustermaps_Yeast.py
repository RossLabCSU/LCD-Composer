
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
organism = 'Scerevisiae'


def main():

    df = {}
    
    for aa in amino_acids:
        df[aa] = df.get(aa, [])
    
    df['Proteins'] = []
    df['Sequences'] = []
    df['LCD Category'] = []
    df['Secondary AA'] = []

    df = get_all_LCDs(df)
    proteins = df.pop('Proteins')
    seqs = df.pop('Sequences')
    proteome_df = pd.DataFrame.from_dict(df)
    
    for aa in amino_acids.replace('W', ''):     # W IS EXCLUDED BECAUSE THERE WERE NO W-RICH PRIMARY LCDs USING THE STANDARD 40% COMPOSITION THRESHOLD
        lcd_df = proteome_df[proteome_df['LCD Category'] == aa]
        lcd_df = lcd_df.drop('LCD Category', axis=1)
        lcd_df = lcd_df.drop('Proteomes', axis=1)
        lcd_df.sort_values(aa, inplace=True, ascending=False)

        master_array = []
        
        for secondary_aa in amino_acids.replace(aa, ''):
            secondary_df = lcd_df[lcd_df['Secondary AA'] == secondary_aa]
            secondary_df = secondary_df.drop('Secondary AA', axis=1)
            secondary_df.sort_values(secondary_aa, inplace=True, ascending=False)
            
            df = pd.DataFrame.to_dict(secondary_df)

            for num in df[secondary_aa]:
                l = []
                for key in df:
                    l.append( df[key][num] )
                master_array.append(l)

        heatmap(master_array, aa)

    
def heatmap(master_array, aa):
    n = len(master_array)
    cg = sns.heatmap(master_array, yticklabels=True, vmin=0, vmax=100, cbar=False)
    plt.xticks([x+0.5 for x in range(20)], labels=list(amino_acids), fontname='Arial', fontsize=18, ha='center')
    plt.xlabel('Amino Acid', fontname='Arial', fontsize=18)
    plt.title(aa + '-rich LCDs', fontname='Arial', fontsize=18)
    plt.ylabel('Individual Domains (n = ' + str(n) + ')', fontname='Arial', fontsize=18)
    ax = plt.gca()
    ax.set_yticklabels([])    

    fig = plt.gcf()

    #====Specific sizing for main body figure=======
    if aa in 'DENQT':
        fig.set_size_inches(5.6, 8)
    #===============================================
    
    #Scaled sizing for supplemental figures===============
    else:
        if n < 25:
            fig.set_size_inches(5.6, n*0.09)
        elif n > 200 and n < 800:
            fig.set_size_inches(5.6, n*0.03)
        elif n > 800:
            fig.set_size_inches(5.6, n*0.01)
        else:
            fig.set_size_inches(5.6, n*0.05)
        

    plt.ylim(n, 0)
    plt.savefig(aa + '-rich LCDs SecondaryAA Heatmap.tiff', bbox_inches='tight', dpi=600)
    plt.close()

    
def get_all_LCDs(df):
    
    df['Proteomes'] = []
    
    for aa in amino_acids:
        h = open(organism + '_' + aa + '_40-CompositionThreshold_20-WindowSize_0p5-LinearDispersionThreshold_LCD-Composer_RESULTS.tsv')
        for i in range(7):
            h.readline()
        header = h.readline()
        
        for line in h:
            id, seq, boundaries, final_comp, final_disp, *remainder = line.rstrip().split('\t')
            junk, uniprot_id, junk = id.split('|')
            
            non_primary_comps = [(seq.count(x) / len(seq) * 100) for x in amino_acids.replace(aa, '')]
            non_primary_aas = amino_acids.replace(aa, '')
            secondary_aa = non_primary_aas[ non_primary_comps.index(max(non_primary_comps)) ]
            
            for res in amino_acids:
                df[res].append( seq.count(res) / len(seq) * 100 )
            df['LCD Category'].append(aa)
            df['Proteins'].append(uniprot_id)
            df['Sequences'].append(seq)
            df['Proteomes'].append(organism)
            
            #==============================================
            #THIS AUTOMATICALLY OMITS ALL LCDs THAT DON'T HAVE A CLEAR SECONDARY AA, EVEN THOUGH SOME LCDs MIGHT HAVE MODERATELY STRONG ENRICHMENT (e.g. > 20%) OF MULTIPLE SUBSIDIARY AAs
            if non_primary_comps.count( max(non_primary_comps) ) < 2:
                df['Secondary AA'].append( secondary_aa )
            else:
                df['Secondary AA'].append( np.nan )
            #==============================================
            
        h.close()
            
    return df
        
if __name__ == '__main__':
    main()