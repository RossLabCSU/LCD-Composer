
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

proteomes = ['UP000002311_Scerevisiae_NoIsoforms', 'UP000001940_6239_Celegans_NoIsoforms', 'UP000000803_7227_Dmelanogaster_NoIsoforms', 'UP000000437_7955_Drerio_NoIsoforms', 'UP000000589_10090_Mmusculus_NoIsoforms', 'UP000005640_9606_Hsapiens_NoIsoforms']
abbrevs = ['Scerevisiae', 'Celegans', 'Dmelanogaster', 'Drerio', 'Mmusculus', 'Hsapiens']
    

amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
aa_names = {'A':'Alanine (A)', 'C':'Cysteine (C)', 'D':'Aspartic Acid (D)', 'E':'Glutamic Acid (E)',
        'F':'Phenylalanine (F)', 'G':'Glycine (G)', 'H':'Histidine (H)', 'I':'Isoleucine (I)', 
        'K':'Lysine (K)', 'L':'Leucine (L)', 'M':'Methionine (M)', 'N':'Asparagine (N)',
        'P':'Proline (P)', 'Q':'Glutamine (Q)', 'R':'Arginine (R)', 'S':'Serine (S)',
        'T':'Threonine (T)', 'V':'Valine (V)', 'W':'Tryptophan (W)','Y':'Tyrosine (Y)'}
        
def main():

    df = {'Proteome':[],
        'GO term':[],
        'Amino Acid':[],
        'GO ID':[],
        'Category':[]
        }
        
    for aa in amino_acids:
        all_go_hits = []
        for i in range(len(proteomes)):
        
            # SKIPS LCD CATEGORY IF NO LCDs OF THAT CATEGORY WERE OBSERVED IN THE ORGANISM
            try:
                h = open(abbrevs[i] + '_' + aa + '_GO_RESULTS_EXCLUDED-ANNOTS.tsv')
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
    plotting(df)
            

def plotting(df):
        
    aa_index = 1
    fig, ax = plt.subplots(5,4,sharex=True, sharey=False, figsize=(8,10))
    big_subplot = fig.add_subplot(111)
    big_subplot.spines['top'].set_color('none')
    big_subplot.spines['bottom'].set_color('none')
    big_subplot.spines['left'].set_color('none')
    big_subplot.spines['right'].set_color('none')
    big_subplot.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
    
    df = pd.DataFrame.from_dict(df)
    
    counts_df = {}
    for aa in amino_acids:
        counts_df[aa] = []
        temp_df = df[df['Amino Acid'] == aa]
        counts = []
        for go_term in list(set(temp_df['GO term'])):
            counts.append( list(temp_df['GO term']).count(go_term) )
            counts_df[aa].append( list(temp_df['GO term']).count(go_term) )

    for aa in amino_acids:
        new_df = {'LCD Class':[],
                'Number of Organisms':[],
                'Count':[],
                'Category':[]
                }
        if len(counts_df[aa]) == 0:
            continue
        xvals = []
        yvals = []
        total_obs = 0
        for x in range(1, 8):
            xvals.append(x)
            yvals.append( counts_df[aa].count(x) )
            total_obs += counts_df[aa].count(x)
            new_df['LCD Class'].append( aa )
            new_df['Number of Organisms'].append( x )
            new_df['Count'].append( counts_df[aa].count(x) )
            new_df['Category'].append( 'Observed' )
            
        #PLOTTING
        ax = plt.subplot(5,4, aa_index)
        sns.barplot(x=xvals, y=yvals, palette=['#1f77b4'])   #PLOT OBSERVED VALUES

        overall_max = ax.patches[0].get_height()
        plt.ylim(0, overall_max + overall_max*0.2)

        plt.ylabel('')
        plt.xlabel('')
        plt.title(aa_names[aa], fontsize=12, fontname='Arial', pad=0.1)
        plt.xlim(-1, 7)
        plt.xticks([x for x in range(8)], labels=[str(x) for x in range(1, 8)], fontname='Arial', fontsize=14)
        
        if aa == 'W':
            plt.yticks([0, 1], labels=['0', '1'], fontname='Arial', fontsize=14)
        elif aa == 'F':
            plt.yticks([0, 2, 4, 6, 8], labels=['0', '2', '4', '6', '8'], fontname='Arial', fontsize=14)
        elif aa == 'M':
            plt.yticks([0, 1, 2], labels=['0', '1', '2'], fontname='Arial', fontsize=14)
        elif aa == 'Y':
            plt.yticks([0, 1, 2, 3, 4], labels=['0', '1', '2', '3', '4'], fontname='Arial', fontsize=14)
        elif aa == 'C':
            plt.yticks([0, 2, 4, 6, 8, 10, 12], labels=['0', '2', '4', '6', '8', '10', '12'], fontname='Arial', fontsize=14)
        else:
            plt.yticks(fontname='Arial', fontsize=14)

        ax = plt.gca()
        if aa_index in range(1,17):
            ax.set_xticklabels([])

        aa_index += 1

    fig.text(-0.025, 0.5, 'Frequency (Number of Enriched GO Terms)', va='center', rotation='vertical', fontname='Arial', fontsize=16)
    fig.text(0.5, -0.02, 'Number of Organisms Sharing Enriched GO Term', ha='center', fontname='Arial', fontsize=16)
    plt.tight_layout(pad=0.2)
    plt.savefig('Fig S10 - ModelEukaryoticOrganisms_EXCLUDED-ANNOTS_ProteinSampling_GOresults.tiff', bbox_inches ='tight', dpi=600)
    plt.close()  
        
    
if __name__ == '__main__':
    main()
    