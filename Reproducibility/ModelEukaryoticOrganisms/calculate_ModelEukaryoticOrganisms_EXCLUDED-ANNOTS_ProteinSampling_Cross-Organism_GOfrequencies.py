
abbrevs = ['Scerevisiae', 'Celegans', 'Dmelanogaster', 'Drerio', 'Mmusculus', 'Hsapiens']
import statistics
import math
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'

aa_names = {'A':'Alanine (A)', 'C':'Cysteine (C)', 'D':'Aspartic Acid (D)', 'E':'Glutamic Acid (E)',
        'F':'Phenylalanine (F)', 'G':'Glycine (G)', 'H':'Histidine (H)', 'I':'Isoleucine (I)', 
        'K':'Lysine (K)', 'L':'Leucine (L)', 'M':'Methionine (M)', 'N':'Asparagine (N)',
        'P':'Proline (P)', 'Q':'Glutamine (Q)', 'R':'Arginine (R)', 'S':'Serine (S)',
        'T':'Threonine (T)', 'V':'Valine (V)', 'W':'Tryptophan (W)','Y':'Tyrosine (Y)'}
        
def main():

    output = open('Combined_GOtermResults_ModelEukaryoticOrganisms_EXCLUDED-ANNOTS_ProteinSampling.tsv', 'w')
    output.write('\t'.join(['# GO','NS','enrichment','name','ratio_in_study','ratio_in_pop','p_uncorrected','depth','study_count','p_bonferroni','p_sidak','p_holm','study_items','Iteration #','LCD Class','Organism']) + '\n')
    
    df = {}
    for organism in abbrevs:
        df[organism] = {}
        for aa in amino_acids:
            df[organism][aa] = {}
            for iter_num in range(1, 1001):
                df[organism][aa][iter_num] = []
            
    for organism in abbrevs:
        for aa in amino_acids:
            df, output = get_RandomGO_Results(organism, df, aa, output)
    output.close()
    
    freqs_out = open('ModelEukaryoticOrganisms_Cross-Organism_GOfrequencies_EXCLUDED-ANNOTS_ProteinSampling.csv', 'w')
    freqs_out.write('LCD Class,Number of Organisms with Shared GO term,Mean Number of GO Terms,Standard Deviation,Standard Error of the Mean,95% CI\n')
    
    for lcd_class in amino_acids:
        num_gos_shared = {1:[],
                    2:[],
                    3:[],
                    4:[],
                    5:[],
                    6:[],
                    7:[]}
        for iter_num in range(1, 1001):
            go_set = set()
            frequencies = []
            for organism in abbrevs:
                for go_term in df[organism][lcd_class][iter_num]:
                    go_set.add(go_term)
            for go_term in go_set:
                freq = 0
                for organism in abbrevs:
                    if go_term in df[organism][lcd_class][iter_num]:
                        freq += 1

                frequencies.append(freq)
                
            for i in range(1,8):
                num_shared = frequencies.count(i)
                num_gos_shared[i].append( num_shared )
                
        for i in range(1,8):
            ave = statistics.mean( num_gos_shared[i] )
            stddev = statistics.stdev( num_gos_shared[i] )
            stderr = stddev / math.sqrt( len(num_gos_shared[i]) )
            ci = 1.96 * stderr
            freqs_out.write(','.join( [str(x) for x in [lcd_class, i, ave, stddev, stderr, ci]] ) + '\n')
            
    freqs_out.close()


def get_RandomGO_Results(organism, df, aa, output):

    h = open(organism + '_RandomlySampledProteins_All_Iterations_GO_RESULTS_EXCLUDED-ANNOTS.tsv')
    header = h.readline()
        
    for line in h:
        items = line.rstrip().split('\t')
        go_id, go_term, sidak_pval, iter_num, lcd_class = items[0], items[3], float(items[10]), int(items[-2]), items[-1]
        go_term = go_term.rstrip()
        if lcd_class != aa:
            continue

        df[organism][lcd_class][iter_num].append(go_term)
        output.write('\t'.join(items) + '\t' + organism + '\n')
    h.close()
        
    return df, output
    

if __name__ == '__main__':
    main()