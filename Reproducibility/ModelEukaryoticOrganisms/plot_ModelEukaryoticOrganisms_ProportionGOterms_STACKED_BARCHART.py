
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn as sns
import pandas as pd
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
aa_names = {'A':'Alanine (A)', 'C':'Cysteine (C)', 'D':'Aspartic Acid (D)', 'E':'Glutamic Acid (E)',
        'F':'Phenylalanine (F)', 'G':'Glycine (G)', 'H':'Histidine (H)', 'I':'Isoleucine (I)', 
        'K':'Lysine (K)', 'L':'Leucine (L)', 'M':'Methionine (M)', 'N':'Asparagine (N)',
        'P':'Proline (P)', 'Q':'Glutamine (Q)', 'R':'Arginine (R)', 'S':'Serine (S)',
        'T':'Threonine (T)', 'V':'Valine (V)', 'W':'Tryptophan (W)','Y':'Tyrosine (Y)'}
        
def main():

    proteomes = ['UP000002311_Scerevisiae_NoIsoforms', 'UP000001940_6239_Celegans_NoIsoforms', 'UP000000803_7227_Dmelanogaster_NoIsoforms', 'UP000000437_7955_Drerio_NoIsoforms', 'UP000186698_8355_Xlaevis_NoIsoforms', 'UP000000589_10090_Mmusculus_NoIsoforms', 'UP000005640_9606_Hsapiens_NoIsoforms']
    abbrevs = ['Scerevisiae', 'Celegans', 'Dmelanogaster', 'Drerio', 'Xlaevis', 'Mmusculus', 'Hsapiens']
    
    df = {'Proteome':[],
        'GO term':[],
        'Amino Acid':[],
        'GO ID':[],
        'Category':[]
        }
        
    for aa in amino_acids:
        all_go_hits = []
        for i in range(len(proteomes)):
            try:
                h = open(abbrevs[i] + '_' + aa + '_GO_RESULTS.tsv')
            except:
                continue
            header = h.readline()
            for line in h:
                items = line.rstrip().split('\t')
                if len(items) == 13:
                    go_id, cat, e_or_p, go_desc, ratio_in_list, ratio_in_pop, uncorr_pval, depth, num_prots, bonf_pval, sidak_pval, holm_pval, assoc_prots = items    #lines with 'p' for e_or_p sometimes do not have any proteins in the last cell
                else:
                    go_id, cat, e_or_p, go_desc, ratio_in_list, ratio_in_pop, uncorr_pval, depth, num_prots, bonf_pval, sidak_pval, holm_pval = items
                
                if e_or_p == 'e' and float(sidak_pval) < 0.05:
                    df['Proteome'].append( abbrevs[i] )
                    df['GO term'].append( go_desc )
                    df['GO ID'].append( go_id )
                    df['Amino Acid'].append( aa )
                    df['Category'].append( cat )
                
            h.close()

    stacked_barchart(df)
            

def stacked_barchart(df):
    df = pd.DataFrame.from_dict(df)
    colors = sns.color_palette("colorblind", 7)
    
    counts_df = {}
    num_gos = []
    for aa in amino_acids:
        counts_df[aa] = []
        temp_df = df[df['Amino Acid'] == aa]
        counts = []
        for go_term in list(set(temp_df['GO term'])):
            counts.append( list(temp_df['GO term']).count(go_term) )
            counts_df[aa].append( list(temp_df['GO term']).count(go_term) )
        num_gos.append(len(temp_df['GO term']))

    num_gos, sorted_aas = zip(*sorted(zip(num_gos, list(amino_acids)), reverse=True))

    x_index = 0
    for aa in sorted_aas:
        if len(counts_df[aa]) == 0:
            continue
        xvals = []
        yvals = []
        for x in range(min(counts_df[aa]), max(counts_df[aa])+1):
            xvals.append(x)
            yvals.append( counts_df[aa].count(x) )
            
        height_sum = 0
        colors_index = 0
        for yval in yvals:
            plt.bar(x=x_index, height=yval/sum(yvals), bottom=height_sum, color=colors[colors_index])
            colors_index += 1
            height_sum += yval/sum(yvals)
        plt.text(x=x_index, y=1.01, s='n = ' + str(sum(yvals)), rotation=90, fontname='Arial', fontsize=14, ha='center', va='bottom')
            
        x_index += 1
        
    plt.xlim(-0.6, 19.6)
    plt.ylim([0, 1.25])
    plt.xticks([x for x in range(len(sorted_aas))], labels=list(sorted_aas), fontname='Arial', fontsize=16)
    plt.yticks([x/5 for x in range(6)], labels=[0.0, 0.2, 0.4, 0.6, 0.8, 1.0], fontname='Arial', fontsize=16)
    plt.xlabel('LCD Class', fontname='Arial', fontsize=18)
    plt.ylabel('Proportion of All\nEnriched GO Terms', fontname='Arial', fontsize=18)
    
    leg_items = [Patch(facecolor=colors[i], label=i+1) for i in range(7)]
    leg = plt.legend(handles=leg_items, prop={'family':'Arial', 'size':14}, loc=2, bbox_to_anchor=(1.0, 1.0), title='Number of Organisms\nwith GO Enriched')
    leg.set_title('Number of Organisms\nwith GO Enriched')
    plt.setp(leg.get_title(), multialignment='center', fontname='Arial', fontsize=14)
    
    fig = plt.gcf()
    fig.set_size_inches(8, 5)
    plt.savefig('Fig S12 - ModelEukaryoticOrganisms - Proportion CrossSpecies GOterm Enrichment.tiff', bbox_inches ='tight', dpi=600)
    plt.close()
    
    
if __name__ == '__main__':
    main()
    