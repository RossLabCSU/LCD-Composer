
import pandas as pd
import numpy as np
from scipy import stats
import math
import pickle
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
bg_freqs = pickle.load(open('Background_X-rich_Window_SecondaryAA_Frequencies.dat', 'rb'))

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
        
    output = open('SecondaryAA_LCD_Statistics_Uncorrected_Pvals.csv', 'w')
    output.write('Fisher\'s exact tests for observed vs. expected frequencies for LCD subclasses within each primary LCD class.\n')
    output.write('Primary LCD Class,Secondary AA,P-value,Odds Ratio,lnOR,# Expected with Secondary AA,Total # of Windows,# Observed with Secondary AA,Total # LCDs with Unambiguous SecondaryAA\n')
    
    residues = amino_acids[:].replace('W', '')
    
    for aa in residues:
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
            
        oddsratios, pvals, output = calc_stats(master_array, proteome, aa, output)

        
def calc_stats(master_array, proteome, primary_aa, output):
    
    secondary_aa_bg_freqs = bg_freqs[primary_aa]
    total = len(master_array)

    contingencies = []
    imputed_contingencies = []
    
    residues = amino_acids.replace(primary_aa, '')
    for aa in residues:
        hits = 0
        for i in range(len(master_array)):
            arr = master_array[i][:amino_acids.index(primary_aa)] + master_array[i][amino_acids.index(primary_aa)+1:]
            for j in range(len(arr)):
                if arr[j] == max(arr) and residues[j] == aa:
                    hits += 1
        
        table = []
        num_expected = secondary_aa_bg_freqs[aa] / secondary_aa_bg_freqs['Total Windows'] * total
        if int(num_expected) >= 1 and hits == 0:
            table.append( [ hits+1 , (total - hits) ] )   #OBSERVED - impute value for biased estimate if # observed is ==0 but # expected is >=1
        else:
            table.append( [ hits , (total - hits) ] )   #OBSERVED
        table.append( [ secondary_aa_bg_freqs[aa], (secondary_aa_bg_freqs['Total Windows'] - secondary_aa_bg_freqs[aa]) ] )     #EXPECTED
        contingencies.append(table)

        oddsratio, pval = stats.fisher_exact(table)
        
        num_obs = hits
        if num_obs == 0 and int(num_expected) < 1:
            output.write(','.join([str(x) for x in [primary_aa, aa, 'nan', 'nan', 'nan', secondary_aa_bg_freqs[aa], secondary_aa_bg_freqs['Total Windows'], num_obs, total]]) + '\n')
        elif  num_obs == 0 and int(num_expected) >= 1:
            output.write(','.join([str(x) for x in [primary_aa, aa, pval, oddsratio, math.log(oddsratio), secondary_aa_bg_freqs[aa], secondary_aa_bg_freqs['Total Windows'], str(num_obs+1) + ' (imputed value)', total]]) + '\n')
        else:
            output.write(','.join([str(x) for x in [primary_aa, aa, pval, oddsratio, math.log(oddsratio), secondary_aa_bg_freqs[aa], secondary_aa_bg_freqs['Total Windows'], num_obs, total]]) + '\n')

    oddsratios = []
    pvals = []
    for contingency in contingencies:
        oddsratio, pval = stats.fisher_exact(contingency)
        oddsratios.append( oddsratio )
        pvals.append( pval )
        
    return oddsratios, pvals, output
        
    
    
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