

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import statistics
from scipy import stats
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
organisms = ['Scerevisiae', 'Celegans', 'Dmelanogaster', 'Drerio', 'Xlaevis', 'Mmusculus', 'Hsapiens']

def main():

    diag_df = {'Position':[],
            'Organism Pair':[],
            'Mean Proportion of Overlap':[],
            'Category':[]}
            
    output = open('TableS8_CrossOrganism_Homotypic_vs_Heterotypic_STATISTICS.tsv', 'w')
    output.write('\t'.join( ['Organism1', 'Organism2', 'On-diagonal Values', 'Off-diagonal Values', 'U-statistic', 'Two-sided P-value'] ) + '\n')
            
    pvals = []
    index = 0
    for organism in organisms:
        df = {aa:[] for aa in amino_acids}
        for aa in amino_acids:
            try:
                h = open(organism + '_' + aa + '_GO_RESULTS.tsv')
            except:
                continue
            df = get_GOterms(aa, organism, df, h)

        for sec_org in [x for x in organisms if x != organism]:
            sec_df = {res:[] for res in amino_acids}
            for res in amino_acids:
                try:
                    h = open(sec_org + '_' + res + '_GO_RESULTS.tsv')
                except:
                    continue
                sec_df = get_GOterms(res, sec_org, sec_df, h)

            transposed_matrix, samp_sizes, other_group_samp_sizes = calc_overlap_matrix(df, sec_df)
            
            ave_diag, stddev_diag, ave_offdiag, stddev_offdiag, diag_df, diag_values, offdiag_values, diagonals_zeros_included, off_diagonals_zeros_included = calc_diagonal_offdiagonal(transposed_matrix, organism, sec_org, samp_sizes, other_group_samp_sizes, diag_df, index)

            ustat, pval = stats.mannwhitneyu(diag_values, offdiag_values, alternative='two-sided')
            pvals.append(pval)

            output.write('\t'.join( [organism, sec_org, ','.join([str(x) for x in diagonals_zeros_included]), ','.join([str(x) for x in off_diagonals_zeros_included]), str(ustat), str(pval)] ) + '\n')

            index += 1
            
    output.close()
    plot_diag_overlap_barplot(diag_df, pvals)
    
                   
def get_GOterms(aa, organism, df, h):

    header = h.readline()
    
    for line in h:
        items = line.rstrip().split('\t')
        e_or_p = items[2]
        go_term = items[3]
        sidak_pval = float(items[10])
        if e_or_p == 'e' and sidak_pval < 0.05 - 0.00000001:
            df[aa].append(go_term)
            
    h.close()
    
    return df
    
def plot_diag_overlap_barplot(diag_df, pvals):
    
    xtick_labels = []
    for organism in organisms:
        for second_organism in [x for x in organisms if x != organism]:
            xtick_labels.append( organism[0] + '. ' + organism[1:] + ' + ' + second_organism[0] + '. ' + second_organism[1:] )
    sns.barplot(x='Position', y='Mean Proportion of Overlap', data=diag_df, hue='Category', ci=None, hue_order=['On Diagonal', 'Off Diagonal'])
    plt.xticks([x for x in range(len(xtick_labels))], labels=xtick_labels, fontname='Arial', fontsize=10, rotation=90)
    plt.yticks(fontname='Arial', fontsize=12)
    
    plt.xlabel('Organism Pair', fontname='Arial', fontsize=14)
    plt.ylabel('Mean Proportion of GO\nTerms that Overlap', fontname='Arial', fontsize=14)
    
    fig = plt.gcf()
    fig.set_size_inches(10, 4.8)
    plt.savefig('Fig S7 - Summarized AcrossOrganism GOterm Overlap_MeanOverlapBarplot.tiff', bbox_inches='tight', dpi=600)
    plt.close()

    
def calc_diagonal_offdiagonal(transposed_matrix, organism, second_organism, samp_sizes, other_group_samp_sizes, diag_df, index):

    diagonals = []
    off_diagonals = []
    diagonals_zeros_included = []
    off_diagonals_zeros_included = []
    # ITERATE OVER ROWS
    for i in range(len(transposed_matrix)):
        # EXCLUDE ROWS WHERE THE GO TERM SAMPLE SIZE IS 0
        if other_group_samp_sizes[i] == 0:  
            for k in range(19):
                off_diagonals_zeros_included.append('*0.0')
            diagonals_zeros_included.append('*0.0')
            continue
        # ITERATE OVER COLUMNS
        for j in range(len(transposed_matrix[i])):
            # EXCLUDE COLUMNS WHERE THE GO TERM SAMPLE SIZE IS 0
            if samp_sizes[j] == 0:
                if i == j:
                    diagonals_zeros_included.append('*0.0')
                else:
                    off_diagonals_zeros_included.append('*0.0')
                continue
                
            if i == j:
                diagonals.append( transposed_matrix[i][j] )
                diag_df['Position'].append(index)
                diag_df['Organism Pair'].append( organism[0] + '. ' + organism[1:] + ' + ' + second_organism[0] + '. ' + second_organism[1:] )
                diag_df['Mean Proportion of Overlap'].append( transposed_matrix[i][j] )
                diag_df['Category'].append( 'On Diagonal' )
                
                diagonals_zeros_included.append( transposed_matrix[i][j] )
            else:
                off_diagonals.append( transposed_matrix[i][j] )
                diag_df['Position'].append(index)
                diag_df['Organism Pair'].append( organism[0] + '. ' + organism[1:] + ' + ' + second_organism[0] + '. ' + second_organism[1:] )
                diag_df['Mean Proportion of Overlap'].append( transposed_matrix[i][j] )
                diag_df['Category'].append( 'Off Diagonal' )
                
                off_diagonals_zeros_included.append( transposed_matrix[i][j] )
                
                
    ave_diag = statistics.mean(diagonals)
    ave_offdiag = statistics.mean(off_diagonals)
            
    stddev_diag = statistics.stdev(diagonals)
    stddev_offdiag = statistics.stdev(off_diagonals)
    
    return ave_diag, stddev_diag, ave_offdiag, stddev_offdiag, diag_df, diagonals, off_diagonals, diagonals_zeros_included, off_diagonals_zeros_included


def calc_overlap_matrix(df, sec_df):

    matrix = []
    samp_sizes = [len(df[aa]) for aa in amino_acids] #Contains the sample sizes for the number of enriched GO terms...used in labels on matrix plot
    other_group_samp_sizes = [len(sec_df[aa]) for aa in amino_acids]
    for aa in amino_acids:
        overlaps = []
        for sec_aa in amino_acids:
            aa_gos = df[aa]
            sec_aa_gos = sec_df[sec_aa]
            
            overlap = 0
            for go_term in aa_gos:
                if go_term in sec_aa_gos:
                    overlap += 1 / len(aa_gos)
            overlaps.append(overlap)

        matrix.append(overlaps)
            
    transposed_matrix = np.array(matrix).T
    
    return transposed_matrix, samp_sizes, other_group_samp_sizes


if __name__ == '__main__':
    main()