
import pickle
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

organism = 'Scerevisiae'
    
def main():
    df = {}

    master_counts = []
    for aa in amino_acids:
        counts = []
        df = {}
        h = open(organism + '_' + aa + '_40-CompositionThreshold_20-WindowSize_0p5-LinearDispersionThreshold_LCD-Composer_RESULTS.tsv')
        for i in range(7):
            h.readline()
        header = h.readline()
        for line in h:
            id, seq, bounds, final_comp, final_disp, *comps = line.rstrip().split('\t')
            junk, uniprot, junk = id.split('|')
            start, end = bounds[1:-1].split('-')
            df[uniprot] = df.get(uniprot, {'LCD Class':[], 'Boundaries':[], 'Sequences':[], 'Lines':[]})
            df[uniprot]['LCD Class'].append(aa)
            df[uniprot]['Boundaries'].append( (int(start), int(end)) )
            df[uniprot]['Sequences'].append(seq)
            df[uniprot]['Lines'].append(','.join( [id, seq, aa, bounds, final_comp, final_disp] + comps) )
        h.close()
        
        for secondary_aa in amino_acids:
            h = open(organism + '_' + secondary_aa + '_40-CompositionThreshold_20-WindowSize_0p5-LinearDispersionThreshold_LCD-Composer_RESULTS.tsv')
            
            for i in range(7):
                h.readline()
            header = h.readline()
            hit_prots = set()
            for line in h:
                id, seq, bounds, final_comp, final_disp, *comps = line.rstrip().split('\t')
                junk, uniprot, junk = id.split('|')
                start, end = bounds[1:-1].split('-')
                if uniprot in df:
                    for existing_bounds in df[uniprot]['Boundaries']:
                        existing_positions = [x for x in range(existing_bounds[0], existing_bounds[1]+1)]
                        current_lcd_positions = [x for x in range(int(start), int(end)+1)]
                        check_overlap = [False if x not in existing_positions else True for x in current_lcd_positions]    #CHECKS FOR ANY OVERLAP BETWEEN TWO DOMAINS
                        if aa == 'A' and secondary_aa == 'Q':
                            print(prot, existing_positions, current_lcd_positions, check_overlap)
                            
                        if True not in check_overlap:
                            hit_prots.add(uniprot)
                        else:
                            continue
                                
            if len(hit_prots) > 0 and aa != secondary_aa:
                output = open(organism + '_' + aa + '_' + secondary_aa + '_MultiLCD_Proteins.txt', 'w')
                for prot in hit_prots:
                    output.write(prot + '\n')
                output.close()
                
            if secondary_aa == aa:
                counts.append( 0 )
            else:
                counts.append( len(hit_prots) )
        master_counts.append( counts )
               
    heatmap(master_counts, organism)

        
def heatmap(master_counts, organism):

    mask = np.zeros_like(master_counts)
    mask[np.triu_indices_from(mask)] = True
    ax = sns.heatmap(master_counts, annot=True, fmt='d', annot_kws={"size": 18, 'family':'Arial'}, mask=mask)

    plt.xticks([x+0.5 for x in range(20)], labels=list(amino_acids), fontname='Arial', fontsize=20)
    plt.yticks([x+0.5 for x in range(20)], labels=list(amino_acids), rotation=0, fontname='Arial', fontsize=20)
    plt.xlabel('Primary LCD Class', fontname='Arial', fontsize=22)
    plt.ylabel('Secondary LCD Class', fontname='Arial', fontsize=22)

    cbar = ax.collections[0].colorbar
    labels = cbar.ax.get_yticklabels()
    cbar.ax.set_yticklabels(labels, fontname='Arial', fontsize=20, rotation=90)
    
    cbar.set_label('Number of Proteins with\nCo-occurring LCDs', fontname='Arial', fontsize=18)

    plt.xlim(0, 19)
    plt.ylim(20, 1)
    
    fig = plt.gcf()
    fig.set_size_inches(11.5, 7)
    plt.savefig('Fig S10B - Scerevisiae LCD Co-occurrence in MultiLCD Proteins_HEATMAP.tiff', bbox_inches='tight', dpi=600)
    plt.close()


if __name__ == '__main__':
    main()
