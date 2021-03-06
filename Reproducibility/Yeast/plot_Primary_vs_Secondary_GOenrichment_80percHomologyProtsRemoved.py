
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
from matplotlib.patches import Patch

def main():

    h = open('Primary_vs_Secondary_GOenrichment_80percHomologyProtsRemoved.tsv')
    
    header = h.readline()
    df = {'GO term':[],
            'LCD Class':[],
            'Effect':[]}
    
    #GETS ALL GO TERMS AND THEIR ASSOCIATED EFFECTS OF SUBCLASSING
    #USED BELOW TO CHECK IF THE GO TERM IS RETAINED IN ANY OF THE ASSOCIATED SUBCLASSES AFTER 80% HOMOLOGY PROTEINS ARE REMOVED
    track_go_terms = {lcd_class:{} for lcd_class in amino_acids}
    for line in h:
        go_term, lcd_class, secondary_aas, effect = line.rstrip().split('\t')
        if effect == 'New':
            continue
        track_go_terms[lcd_class][go_term] = track_go_terms.get(go_term, [])
        track_go_terms[lcd_class][go_term] = track_go_terms.get(go_term, [])
        track_go_terms[lcd_class][go_term].append(effect)
    h.close()

    
    for lcd_class in track_go_terms:    
        for go_term in track_go_terms[lcd_class]:
            df['GO term'].append(go_term)
            df['LCD Class'].append(lcd_class)
            
            if 'Retained' in track_go_terms[lcd_class][go_term]:
                df['Effect'].append('Retained')
            else:
                df['Effect'].append('Lost')
 
    df = pd.DataFrame.from_dict(df)
    stacked_barchart(df)
        
    
def stacked_barchart(df):

    colors = ['#0485d1', '0.7']

    all_lcds_counts = [0, 0]
    counts_df = {}
    num_gos = []
    for aa in amino_acids:
        counts_df[aa] = []
        temp_df = df[df['LCD Class'] == aa]
        counts = []
        ind = 0
        for effect in ['Retained', 'Lost']:
            counts.append( list(temp_df['Effect']).count( effect ) )
            counts_df[aa].append( list(temp_df['Effect']).count( effect ) )
            all_lcds_counts[ind] += list(temp_df['Effect']).count( effect )
            ind += 1
        num_gos.append(len(temp_df['GO term']))
        
    num_gos, sorted_aas = zip(*sorted(zip(num_gos, list(amino_acids)), reverse=True))

    x_index = 0
    for aa in sorted_aas:
        if len(counts_df[aa]) == 0:
            continue

        yvals = list(counts_df[aa])
        if sum(yvals) == 0:
            continue
            
        height_sum = 0
        colors_index = 0
        for yval in yvals:
            plt.bar(x=x_index, height=yval/sum(yvals), bottom=height_sum, color=colors[colors_index])
            colors_index += 1
            height_sum += yval/sum(yvals)
        plt.text(x=x_index, y=1.01, s='n = ' + str(sum(yvals)), rotation=90, fontname='Arial', fontsize=14, ha='center', va='bottom')
            
        x_index += 1
        
    ind = 0
    height_sum = 0
    colors_index = 0
    for effect in ['Retained', 'Lost']:
        plt.bar(x=x_index, height=all_lcds_counts[ind]/sum(all_lcds_counts), bottom=height_sum, color=colors[colors_index])
        height_sum += all_lcds_counts[ind] / sum(all_lcds_counts)
        ind += 1
        colors_index += 1
    plt.text(x=x_index, y=1.01, s='n = ' + str(sum(all_lcds_counts)), rotation=90, fontname='Arial', fontsize=14, ha='center', va='bottom')
    
        
    plt.xlim(-0.6, 14.6)
    plt.ylim([0, 1.25])
    plt.xticks([x for x in range(len(sorted_aas[:-6])+1)], labels=list(sorted_aas)[:-6] + ['All LCDs'], fontname='Arial', fontsize=16)
    plt.gca().get_xticklabels()[-1].set_rotation(90)
    plt.yticks([x/5 for x in range(6)], labels=[0.0, 0.2, 0.4, 0.6, 0.8, 1.0], fontname='Arial', fontsize=16)
    plt.xlabel('Primary LCD Class', fontname='Arial', fontsize=18)
    plt.ylabel('Proportion of Enriched\nGO Terms', fontname='Arial', fontsize=18)

    leg_items = [Patch(facecolor=colors[i], label=['Retained', 'Lost'][i]) for i in range(2)]
    leg = plt.legend(handles=leg_items, prop={'family':'Arial', 'size':14}, loc=2, bbox_to_anchor=(1.0, 1.0), title='Effect of Subclassing\non GO term')
    leg.set_title('Effect of Removing\n>80% Homology\nProteins')
    plt.setp(leg.get_title(), multialignment='center', fontname='Arial', fontsize=14)
    
    fig = plt.gcf()
    fig.set_size_inches(5, 5)
    plt.savefig('Fig S15B - Primary_vs_Secondary_GOenrichment_Comparison_80perHomologyProtsRemoved.tiff', bbox_inches ='tight', dpi=600)
    plt.close()


if __name__ == '__main__':
    main()