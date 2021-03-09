
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import random
import pickle
import datetime
from scipy import stats

proteomes = ['UP000002311_Scerevisiae_NoIsoforms', 'UP000001940_6239_Celegans_NoIsoforms', 'UP000000803_7227_Dmelanogaster_NoIsoforms', 'UP000000437_7955_Drerio_NoIsoforms', 'UP000186698_8355_Xlaevis_NoIsoforms', 'UP000000589_10090_Mmusculus_NoIsoforms', 'UP000005640_9606_Hsapiens_NoIsoforms']
abbrevs = ['Scerevisiae', 'Celegans', 'Dmelanogaster', 'Drerio', 'Xlaevis', 'Mmusculus', 'Hsapiens'] 

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
        
    rand_samp_counts = random_go_sampling()   #USE ONLY IF YOU WANT TO RE-RUN THE TIME-CONSUMING SAMPLING. OTHERWISE USE THE LINE BELOW WHICH LOADS THE RESULTS OF THE FIRST 100k SAMPLING
    # rand_samp_counts = pickle.load(open('RandomGOtermSampling_CrossOrganismFrequencyEstimationResults_100000_Iterations_ModelEukaryoticOrganisms.dat', 'rb'))     #USE THIS LINE INSTEAD IF YOU HAVE ALREADY PERFORMED THE 100k SAMPLING...THE SAMPLING RESULTS SHOULD BE SAVED IN THIS DATA FILE.

    stats_output = open('TableS9_Statistics_RandomGOtermSampling_ModelEukaryoticOrganisms.csv', 'w')
    stats_output.write('LCD Class,# of Organisms with Shared GO Term,# Observed,Total Observed,# from Random Sampling,Total from Random Sampling,Odds Ratio,p-value,Bonferroni-corrected p-value\n')

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
            
        total_sampling = sum( [rand_samp_counts[aa][num_organisms] for num_organisms in rand_samp_counts[aa]] )
        yvals_randomsampling = []
        for num_organisms in range(1, 8):
            count = rand_samp_counts[aa][num_organisms]
            yvals_randomsampling.append( count )
            new_df['LCD Class'].append( aa )
            new_df['Number of Organisms'].append( num_organisms )
            new_df['Count'].append( count/total_sampling * total_obs )
            new_df['Category'].append( 'Bootstrapping' )

        #CONTINGENCY AND FISHERS EXACT
        corrected_pvals = []
        colors = []
        for i in range(1,8):
            contingency = []
            hits_obs = yvals[i-1]
            hits_sampling = yvals_randomsampling[i-1]
            
            contingency = [ [hits_obs, total_obs-hits_obs]  ,  [hits_sampling, total_sampling-hits_sampling] ]
            
            oddsratio, pval = stats.fisher_exact(contingency)
            corr_pval = min(1, pval*7)
            corrected_pvals.append( corr_pval )
            stats_output.write(','.join( [str(x) for x in (aa, i, hits_obs, total_obs, hits_sampling, total_sampling, oddsratio, pval, corr_pval)] ) + '\n')
            
            if (hits_obs/total_obs) > (hits_sampling/total_sampling) and corr_pval < 0.05:
                colors.append('#0485d1')
            elif (hits_obs/total_obs) < (hits_sampling/total_sampling) and corr_pval < 0.05:
                colors.append('#ff474c')
            else:
                colors.append('#7f7f7f')
            
        #PLOTTING
        ax = plt.subplot(5,4, aa_index)
        sns.barplot(x=xvals, y=yvals, palette=colors)
        ind = 0
        overall_max = ax.patches[0].get_height()
        for i in range(0, 7):
            max_height = ax.patches[i].get_height()
            if corrected_pvals[ind] < 0.05 and corrected_pvals[ind] >= 0.01:
                ax.text(ax.patches[i].get_x()+0.5, max_height + (overall_max*0.05), '*', ha="center", rotation=90) 
            elif corrected_pvals[ind] < 0.01 and corrected_pvals[ind] >= 0.001:
                ax.text(ax.patches[i].get_x()+0.5, max_height + (overall_max*0.05), '**', ha="center", rotation=90) 
            elif corrected_pvals[ind] < 0.001:
                ax.text(ax.patches[i].get_x()+0.5, max_height + (overall_max*0.05), '***', ha="center", rotation=90) 
            ind += 1
                
        plt.ylim(0, overall_max + overall_max*0.2)
        plt.ylabel('')
        plt.xlabel('')

        plt.title(aa_names[aa], fontsize=12, fontname='Arial', pad=0.1)
        plt.xlim(-1, 7)
        plt.xticks([x for x in range(8)], labels=[str(x) for x in range(1, 8)], fontname='Arial', fontsize=14)
        
        if aa == 'W':
            plt.yticks([0, 1], labels=['0', '1'], fontname='Arial', fontsize=14)
        elif aa == 'F':
            plt.yticks([0, 2, 4, 6, 8, 10], labels=['0', '2', '4', '6', '8', '10'], fontname='Arial', fontsize=14)
        elif aa == 'M':
            plt.yticks([0, 1, 2], labels=['0', '1', '2'], fontname='Arial', fontsize=14)
        else:
            plt.yticks(fontname='Arial', fontsize=14)
   
        ax = plt.gca()
        if aa_index in range(1,17):
            ax.set_xticklabels([])
       
        aa_index += 1

    stats_output.close()
    fig.text(-0.025, 0.5, 'Frequency (Number of Enriched GO Terms)', va='center', rotation='vertical', fontname='Arial', fontsize=16)
    fig.text(0.5, -0.02, 'Number of Organisms Sharing Enriched GO Term', ha='center', fontname='Arial', fontsize=16)
    plt.tight_layout(pad=0.2)
    plt.savefig('Fig S11 - CrossSpecies GO term enrichment__RandomGOtermSampling_ModelEukaryoticOrganisms.tiff', bbox_inches ='tight', dpi=600)
    plt.close()  
    
def random_go_sampling():

    org_goterms = pickle.load(open('GOterms_dict_with_ActualAssociatedProtSet_as_Keys.dat', 'rb'))
    thresholds_df = pickle.load(open('OrganismSpecific_and_LCDspecific_NumAssocProteins_Thresholds_for_GOsignificance.dat', 'rb'))
    goterm_counts = get_goterm_counts()
    iterations = 100000
    cross_spec_counts = {}
    for aa in amino_acids:
        cross_spec_counts[aa] = {1:0,
                        2:0,
                        3:0,
                        4:0,
                        5:0,
                        6:0,
                        7:0
                        }
        print(aa, 'Start time:', str(datetime.datetime.now()))
        for i in range(iterations):
            # print('Iteration Number:', i, str(datetime.datetime.now()))
            organism_samples = {}
            go_terms = []
            for organism in abbrevs:
                if aa not in thresholds_df[organism]:
                    continue
                organism_samples[organism] = []
                organism_gos = list( org_goterms[organism].items() )
                count = goterm_counts[aa][organism]
                num_assoc_prots_threshold = thresholds_df[organism][aa]

                # NEW SAMPLING TO GET GO TERMS WITH A NUMBER OF ASSOCIATED PROTEINS ABOVE EACH ORGANISM/LCD CLASS-SPECIFIC THRESHOLD TO ACHIEVE SIGNIFICANCE
                go_samp = []
                while len(go_samp) < count:
                    go_id, prot_set = random.choice( organism_gos )
                    if go_id not in go_samp and len(prot_set) >= num_assoc_prots_threshold:
                        go_samp.append( go_id )

                organism_samples[organism].append( go_samp )
                go_terms += go_samp
                
            for go_term in set(go_terms):
                count = go_terms.count(go_term)
                cross_spec_counts[aa][count] += 1
                
    pf = open('RandomGOtermSampling_CrossOrganismFrequencyEstimationResults_' + str(iterations) + '_Iterations_ModelEukaryoticOrganisms.dat', 'wb')
    pickle.dump(cross_spec_counts, pf)
    pf.close()

    return cross_spec_counts
    
                
def get_goterm_counts():
    df = {}
    for aa in amino_acids:
        df[aa] = {}
        for abbrev in abbrevs:
            try:
                h = open(abbrev + '_' + aa + '_GO_RESULTS.tsv')
            except:
                df[aa][abbrev] = 0
                continue
            header = h.readline()
            count = 0
            for line in h:
                items = line.rstrip().split('\t')
                enr_or_pur = items[2]
                p_sidak = float(items[10])
                if enr_or_pur == 'e' and p_sidak < 0.05 + 0.0000000001:
                    count += 1
            h.close()
            
            df[aa][abbrev] = count
            
    return df
    
    
if __name__ == '__main__':
    main()
    