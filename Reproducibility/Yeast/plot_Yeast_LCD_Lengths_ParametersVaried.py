
import matplotlib.pyplot as plt
import seaborn as sns
import statistics
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
import pandas as pd
aa_names = {'A':'Alanine (A)', 'C':'Cysteine (C)', 'D':'Aspartic Acid (D)', 'E':'Glutamic Acid (E)',
        'F':'Phenylalanine (F)', 'G':'Glycine (G)', 'H':'Histidine (H)', 'I':'Isoleucine (I)', 
        'K':'Lysine (K)', 'L':'Leucine (L)', 'M':'Methionine (M)', 'N':'Asparagine (N)',
        'P':'Proline (P)', 'Q':'Glutamine (Q)', 'R':'Arginine (R)', 'S':'Serine (S)',
        'T':'Threonine (T)', 'V':'Valine (V)', 'W':'Tryptophan (W)','Y':'Tyrosine (Y)'}

def main():

    for category in ['Window Size', 'Dispersion Threshold', 'Composition Threshold']:
        domain_lengths = {'Composition Threshold':[],
                        'Window Size':[],
                        'Dispersion Threshold':[],
                        'Length':[],
                        'Amino Acid':[],
                        'Varied Parameter':[]}
                        
        domain_lengths = get_my_data(domain_lengths, category)

        if category == 'Window Size':
            len_df = pd.DataFrame.from_dict(domain_lengths)
            len_df = len_df.loc[(len_df['Dispersion Threshold']=='0.5') & (len_df['Composition Threshold']==40)]
        elif category == 'Dispersion Threshold':
            len_df = pd.DataFrame.from_dict(domain_lengths)
            len_df = len_df.loc[(len_df['Window Size']==20) & (len_df['Composition Threshold']==40)]
        else:
            len_df = pd.DataFrame.from_dict(domain_lengths)
            len_df = len_df.loc[(len_df['Window Size']==20) & (len_df['Dispersion Threshold']=='0.5')]
        
        huesplit_barcharts(len_df, category)


def get_my_data(domain_lengths, category):
                
    for aa in amino_acids:
        if category == 'Window Size':
            files = ['Scerevisiae_'+aa+'_40-CompositionThreshold_'+str(win_size)+'-WindowSize_0p5-LinearDispersionThreshold_LCD-Composer_RESULTS.tsv' for win_size in range(10, 110, 10)]
            comp_thresholds = [40 for i in range(len(files))]
            disp_thresholds = ['0p5' for i in range(len(files))]
            win_sizes = [x for x in range(10, 110, 10)]
        elif category == 'Dispersion Threshold':
            files = ['Scerevisiae_'+aa+'_40-CompositionThreshold_20-WindowSize_'+disp_threshold+'-LinearDispersionThreshold_LCD-Composer_RESULTS.tsv' for disp_threshold in ['0p0', '0p1', '0p2', '0p3', '0p4', '0p5', '0p6', '0p7', '0p8', '0p9', '1p0']]
            comp_thresholds = [40 for i in range(len(files))]
            disp_thresholds = [x for x in ['0p0', '0p1', '0p2', '0p3', '0p4', '0p5', '0p6', '0p7', '0p8', '0p9', '1p0']]
            win_sizes = [20 for i in range(len(files))]
        else:
            files = ['Scerevisiae_'+aa+'_'+str(comp_threshold)+'-CompositionThreshold_20-WindowSize_0p5-LinearDispersionThreshold_LCD-Composer_RESULTS.tsv' for comp_threshold in range(10, 110, 10)]
            comp_thresholds = [x for x in range(10, 110, 10)]
            disp_thresholds = ['0p5' for i in range(len(files))]
            win_sizes = [20 for i in range(len(files))]
                

        for i in range(len(files)):
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

            domain_lengths['Composition Threshold'].append( comp_thresholds[i]  )
            domain_lengths['Window Size'].append( win_sizes[i]  )
            domain_lengths['Dispersion Threshold'].append( disp_thresholds[i].replace('p', '.')  )
            domain_lengths['Amino Acid'].append( aa )
            domain_lengths['Varied Parameter'].append( category )
            if len(lengths) > 0:
                domain_lengths['Length'].append( statistics.mean(lengths) )
            else:
                domain_lengths['Length'].append( 0 )
                
            h.close()
                
    return domain_lengths
    
    
def huesplit_barcharts(domain_lengths, category):
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf', '0.8']

    sns.barplot(x=domain_lengths['Amino Acid'], y=domain_lengths['Length'], hue=domain_lengths[category], palette=colors)
    
    plt.xticks(fontname='Arial', fontsize=18)
    plt.yticks(fontname='Arial', fontsize=18)
    plt.xlabel('Amino Acid', fontname='Arial', fontsize=20)
    plt.ylabel('Average Domain Length', fontname='Arial', fontsize=20)
    
    leg = plt.legend(prop={'family':'Arial', 'size':14}, loc=2, bbox_to_anchor=(1.0, 1.02))
    leg.set_title(category.replace(' ', '\n'), prop={'family':'Arial', 'size':14})
    plt.tight_layout()
    fig = plt.gcf()
    fig.set_size_inches(12, 6)
    if category == 'Window Size':
        plt.savefig('Fig S5A - AveDomainLengths_' + category + ' Varied.tiff', bbox_inches='tight', dpi=600)
    elif category == 'Dispersion Threshold':
        plt.savefig('Fig S5B - AveDomainLengths_' + category + ' Varied.tiff', bbox_inches='tight', dpi=600)
    else:
        plt.savefig('Fig S5C - AveDomainLengths_' + category + ' Varied.tiff', bbox_inches='tight', dpi=600)
    plt.close()
    

if __name__ == '__main__':
    main()