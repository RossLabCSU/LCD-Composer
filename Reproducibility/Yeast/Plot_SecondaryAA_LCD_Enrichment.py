
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from scipy import stats
import math
import pickle
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
bg_freqs = pickle.load(open('Background_X-rich_Window_SecondaryAA_Frequencies.dat', 'rb'))
  
aa_titles = {'A':'A-rich', 'C':'Cysteine (C)', 'D':'Aspartic Acid (D)', 'E':'Glutamic Acid (E)',
        'F':'Phenylalanine (F)', 'G':'Glycine (G)', 'H':'Histidine (H)', 'I':'Isoleucine (I)', 
        'K':'Lysine (K)', 'L':'Leucine (L)', 'M':'Methionine (M)', 'N':'Asparagine (N)',
        'P':'Proline (P)', 'Q':'Glutamine (Q)', 'R':'Arginine (R)', 'S':'Serine (S)',
        'T':'Threonine (T)', 'V':'Valine (V)', 'W':'Tryptophan (W)','Y':'Tyrosine (Y)'}

def main():
    
    proteome = 'Scerevisiae'

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

    residues = amino_acids[:]
    
    aa_index = 1
    fig, axes = plt.subplots(5,4,sharex=False, sharey=False, figsize=(8,10))
    big_subplot = fig.add_subplot(111)
    big_subplot.spines['top'].set_color('none')
    big_subplot.spines['bottom'].set_color('none')
    big_subplot.spines['left'].set_color('none')
    big_subplot.spines['right'].set_color('none')
    big_subplot.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
    
    for aa in residues:
        if aa == 'W':
            continue

        lcd_df = proteome_df[proteome_df['LCD Category'] == aa]
        lcd_df = lcd_df.drop('LCD Category', axis=1)
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

        oddsratios, pvals, green_positions, orange_positions, num_lcds = calc_stats(master_array, proteome, aa)

        bonf_pvals = [x*(19-len(green_positions)) for x in pvals]  #ONLY 19 POSSIBLE SUBCATEGORIES SINCE PRIMARY AA CANNOT BE USED AS A SUBCLASS
        lnORs = [math.log(x) if x != 0 else 0 for x in oddsratios]   #USE THIS LINE IF AA COLORING ISN'T NECESSARY

        ax = plt.subplot(5,4, aa_index)
        
        sns.barplot(x=[x for x in range(19)], y=lnORs, color='#1f77b4')

        cleaned_lnORs = [x for x in lnORs if not math.isnan(x) and not math.isinf(x)]

        ymin = min(lnORs) - 1
        ymax = max(lnORs) + 1.2
        
        for i in range(19):
            if lnORs[i] < 0:
                offset = -0.3
            else:
                offset = +0.15
            if bonf_pvals[i] < 0.001:
                if lnORs[i] > 0:
                    plt.text(i+0.4, lnORs[i] + (ymax-ymin)*0.09, '***', va='center', ha='center', rotation=90)
                else:
                    plt.text(i+0.4, lnORs[i] - (ymax-ymin)*0.085, '***', va='center', ha='center', rotation=90)
            elif bonf_pvals[i] < 0.01:
                if lnORs[i] > 0:
                    plt.text(i+0.4, lnORs[i] + (ymax-ymin)*0.07, '**', va='center', ha='center', rotation=90)
                else:
                    plt.text(i+0.4, lnORs[i] - (ymax-ymin)*0.055, '**', va='center', ha='center', rotation=90)
            elif bonf_pvals[i] < 0.05 and lnORs[i] != 0:
                if lnORs[i] > 0:
                    plt.text(i+0.4, lnORs[i] + (ymax-ymin)*0.04, '*', va='center', ha='center', rotation=90)
                else:
                    plt.text(i+0.4, lnORs[i] - (ymax-ymin)*0.03, '*', va='center', ha='center', rotation=90)

        plt.title(aa + '-rich LCDs (n=' + str(num_lcds) + ')', fontsize=12, fontname='Arial', pad=0)
        plt.xticks([x for x in range(19)], labels=list(amino_acids.replace(aa, '')), fontname='Arial', fontsize=7.5)
        
        for x in green_positions:
            ax.get_xticklabels()[x].set_color('#029386')
        for x in orange_positions:
            ax.get_xticklabels()[x].set_color('#ff7f0e')
            
        if aa in 'DN':
            plt.ylim(min(cleaned_lnORs) - 1.35, max(cleaned_lnORs) + 1.35)
        else:
            plt.ylim(min(cleaned_lnORs) - 1, max(cleaned_lnORs) + 1.2)
        plt.yticks(fontname='Arial', fontsize=11)
  
        plt.plot([-1, 20], [0, 0], color='0.2', linewidth=1)
        plt.xlim(-0.75, 19.75)
        plt.tight_layout(w_pad=0, h_pad=0.7)
        
        aa_index += 1

    fig.delaxes(ax = axes[-1,-1])
    fig.text(-0.01, 0.5, 'lnOR', va='center', rotation='vertical', fontname='Arial', fontsize=16)
    fig.text(0.5, -0.01, 'Secondary AA', ha='center', fontname='Arial', fontsize=16)
    plt.savefig('Fig 8 - Scerevisiae_SecondaryAA_Enrichment_Plots.tiff', bbox_inches='tight', dpi=600)
    plt.close()

        
def calc_stats(master_array, proteome, primary_aa):
        
    secondary_aa_bg_freqs = bg_freqs[primary_aa]
    
    total = len(master_array)

    contingencies = []
    imputed_contingencies = []
    
    green_positions = []
    orange_positions = []
    
    index = 0
    residues = amino_acids.replace(primary_aa, '')
    for aa in residues:
        hits = 0
        for i in range(len(master_array)):
            arr = master_array[i][:amino_acids.index(primary_aa)] + master_array[i][amino_acids.index(primary_aa)+1:]
            for j in range(len(arr)):
                if arr[j] == max(arr) and residues[j] == aa:
                    hits += 1
        
        table = []
        table.append( [ hits , (total - hits) ] )   #OBSERVED
        table.append( [ secondary_aa_bg_freqs[aa], (secondary_aa_bg_freqs['Total Windows'] - secondary_aa_bg_freqs[aa]) ] )     #EXPECTED
        contingencies.append(table)
        
        
        num_expected = secondary_aa_bg_freqs[aa] / secondary_aa_bg_freqs['Total Windows'] * total
        if int(num_expected) >= 1 and hits == 0:
            observed = [ hits+1 , (total - hits) ]
            expected = [ secondary_aa_bg_freqs[aa], (secondary_aa_bg_freqs['Total Windows'] - secondary_aa_bg_freqs[aa]) ]     #EXPECTED
            table = [observed, expected]
            orange_positions.append(index)
        elif num_expected < 1 and hits == 0:
            green_positions.append(index)

        imputed_contingencies.append(table)
        index += 1

    oddsratios = []
    pvals = []
    for i in range(len(contingencies)):
        imputed_oddsratio, imputed_pval = stats.fisher_exact(imputed_contingencies[i])
        oddsratio, pval = stats.fisher_exact(contingencies[i])
        oddsratios.append( imputed_oddsratio )
        pvals.append( pval )
        
    return oddsratios, pvals, green_positions, orange_positions, total
        
    
def get_all_LCDs(df):
  
    for aa in amino_acids:
        h = open('Scerevisiae_' + aa + '_40-CompositionThreshold_20-WindowSize_0p5-LinearDispersionThreshold_LCD-Composer_RESULTS.tsv')
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
            
            #==============================================
            #THIS AUTOMATICALLY OMITS ALL LCDs THAT DON'T HAVE A CLEAR SECONDARY AA, EVEN THOUGH SOME LCDs MIGHT HAVE MODERATELY STRONG ENRICHMENT (e.g. > 20%) OF MULTIPLE SUBSIDIARY AAs
            if non_primary_comps.count( max(non_primary_comps) ) < 2:
                df['Secondary AA'].append( secondary_aa )
            else:
                df['Secondary AA'].append( np.nan )
            #==============================================
            
    return df
        
if __name__ == '__main__':
    main()