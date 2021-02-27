
import matplotlib.pyplot as plt
import seaborn as sns
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
import pandas as pd
import numpy as np
from Bio import SeqIO
aa_names = {'A':'Alanine (A)', 'C':'Cysteine (C)', 'D':'Aspartic Acid (D)', 'E':'Glutamic Acid (E)',
        'F':'Phenylalanine (F)', 'G':'Glycine (G)', 'H':'Histidine (H)', 'I':'Isoleucine (I)', 
        'K':'Lysine (K)', 'L':'Leucine (L)', 'M':'Methionine (M)', 'N':'Asparagine (N)',
        'P':'Proline (P)', 'Q':'Glutamine (Q)', 'R':'Arginine (R)', 'S':'Serine (S)',
        'T':'Threonine (T)', 'V':'Valine (V)', 'W':'Tryptophan (W)','Y':'Tyrosine (Y)'}

def main():

    for category in ['Window Size', 'Dispersion Threshold', 'Composition Threshold']:
        total_domains_wideform = {} 
        total_domains_wideform = get_my_data(category)

        total_domains_wideform = pd.DataFrame.from_dict(total_domains_wideform)

        heatmap_total_nums(total_domains_wideform, category)
    
    
def get_my_data(category):

    proteome_len = get_proteome_length()
    
    totals_wideform = {'Amino Acid':[]}
    for aa in amino_acids:
        if category == 'Window Size':
            files = ['Scerevisiae_'+aa+'_40-CompositionThreshold_'+str(win_size)+'-WindowSize_0p5-LinearDispersionThreshold_LCD-Composer_RESULTS.tsv' for win_size in range(10, 110, 10)]
            varied_params = [x for x in range(10, 110, 10)]
        elif category == 'Dispersion Threshold':
            files = ['Scerevisiae_'+aa+'_40-CompositionThreshold_20-WindowSize_'+disp_threshold+'-LinearDispersionThreshold_LCD-Composer_RESULTS.tsv' for disp_threshold in ['0p0', '0p1', '0p2', '0p3', '0p4', '0p5', '0p6', '0p7', '0p8', '0p9', '1p0']]
            varied_params = [x for x in ['0.0', '0.1', '0.2', '0.3', '0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1.0']]
        else:
            files = ['Scerevisiae_'+aa+'_'+str(comp_threshold)+'-CompositionThreshold_20-WindowSize_0p5-LinearDispersionThreshold_LCD-Composer_RESULTS.tsv' for comp_threshold in range(10, 110, 10)]
            varied_params = [x for x in range(10, 110, 10)]
            
        totals_wideform['Amino Acid'].append(aa)
        for i in range(len(files)):
            totals_wideform[str(varied_params[i])] = totals_wideform.get(str(varied_params[i]), [])
            h = open(files[i])
            for j in range(7):
                h.readline()
            header = h.readline()
            
            counts = 0
            lengths = []
            for line in h:
                id, seq, boundaries, final_comp, final_stddev, *remainder = line.rstrip().split('\t')
                counts += 1
                lengths.append(len(seq))
                
            h.close()

            if category in ['Window Size', 'Dispersion Threshold']:
                totals_wideform[str(varied_params[i])].append( counts )
            else:
                totals_wideform[str(varied_params[i])].append( round(sum(lengths) / proteome_len * 100, 2) )

        h.close()

    return totals_wideform
    
    
def heatmap_total_nums(df, category):
    
    df.set_index('Amino Acid', inplace=True)

    if category == 'Window Size':
        df.sort_values('20', inplace=True, ascending=False)
    elif category == 'Composition Threshold':
        df.sort_values('10', inplace=True, ascending=False)
    else:
        df.sort_values('0.0', inplace=True, ascending=False)
        
    colors = sns.color_palette("coolwarm", 1000)
    
    log_df = np.log10(df.replace(0, np.nan))
    if category in ['Window Size', 'Dispersion Threshold']:
        log_df.replace(np.nan, 0, inplace=True)
    else:
        log_df.replace(np.nan, -2.1, inplace=True)  # -2.1 REPRESENTS THE MINIMUM OBSERVED log(fraction proteome in LCDs) FOR HEATMAP COLOR SCALE

    ax = sns.heatmap(log_df, cmap=colors, annot=df, fmt='#', annot_kws={'family':'Arial', 'size':12}, cbar=False)

    plt.xticks(fontname='Arial', fontsize=16)
    plt.yticks(fontname='Arial', fontsize=16, rotation=0)
    plt.xlabel('Linear Dispersion Threshold', fontname='Arial', fontsize=18)
    plt.ylabel('Amino Acid', fontname='Arial', fontsize=18)

    plt.ylim(20, 0)
    fig = plt.gcf()
    fig.set_size_inches(5.5, 8)
    if category == 'Window Size':
        plt.savefig('Fig S4A - TotalNumberOfLCDs_' + category + ' Varied.tiff', bbox_inches='tight', dpi=600)
    elif category == 'Dispersion Threshold':
        plt.savefig('Fig S4B - TotalNumberOfLCDs_' + category + ' Varied.tiff', bbox_inches='tight', dpi=600)
    else:
        plt.savefig('Fig S4C - TotalNumberOfLCDs_' + category + ' Varied.tiff', bbox_inches='tight', dpi=600)
    plt.close()
    
    
def get_proteome_length():
    h = open('UP000002311_Scerevisiae_NoIsoforms.fasta')
    total_len = 0
    for seq_record in SeqIO.parse(h, 'fasta'):
        seq = str(seq_record.seq)
        if seq[-1] == '*':
            seq = seq[:-1]
            
        total_len += len(seq)
        
    return total_len
    

if __name__ == '__main__':
    main()