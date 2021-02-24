

amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
import numpy as np
from Bio import SeqIO
import pickle
taxonid_to_organism = pickle.load(open('TaxonID_to_Organism (2).dat', 'rb'))    #NEED TO CHANGE BACK AFTER TESTING

def main():

    proteomes = get_proteomes()
    df = {'Proteome':[],
        'LCD Class':[],
        'Minimum':[],
        '25th Percentile':[],
        '50th Percentile':[],
        '75th Percentile':[],
        'Maximum':[]}
    output = open('Viruses_LCDcontent_Statistics_TESTING.tsv', 'w') #NEED TO CHANGE BACK
    output.write('All statistics (except Total Percentage of Proteome in LCD) are per-protein LCD percentages for each LCD class (i.e. the amount of the sequence that is covered by the given LCD type)\n')
    output.write('\t'.join( ['Proteome', 'Organism Name', 'LCD Class', 'Minimum', '25th Percentile', '50th Percentile (Median)', '75th Percentile', 'Maximum', 'Average', 'Total Percentage of Proteome in LCD']) + '\n')
    
    for proteome in proteomes:
        #ONLY SKIPS _additional FILES FOR THE INITIAL SEARCH...EACH _additional FILE IS CHECKED FOR AND INCORPORATED INTO THE DATA FOR THE CORRESPONDING ORIGINAL PROTEOME IN THE CHECKS BELOW
        if '_additional' in proteome:
            continue
        
        prots, seqs_df = get_prots(proteome)
        prot_df = {prot:{aa:0 for aa in amino_acids} for prot in prots}
        lcd_content_df = {aa:[] for aa in amino_acids}
        concatenated_lcd_length = {aa:0 for aa in amino_acids}
        
        total_proteome_length = sum( [len(seqs_df[prot]) for prot in seqs_df] )
                
        prot_df, lcd_content_df, concatenated_lcd_lengths = calculate_perc_lcd(proteome, prot_df, lcd_content_df, seqs_df, concatenated_lcd_length)
        
        if proteome + '_additional' in proteomes:
            prots, seqs_df = get_prots(proteome + '_additional')
            for prot in prots:
                prot_df[prot] = {aa:0 for aa in amino_acids}
            prot_df, lcd_content_df, concatenated_lcd_lengths = calculate_perc_lcd(proteome + '_additional', prot_df, lcd_content_df, seqs_df, concatenated_lcd_length)
            
        uniprot, taxonid = proteome.split('_')
        organism = taxonid_to_organism[ taxonid ]
                    
        for aa in amino_acids:
            arr = np.array(lcd_content_df[aa])
            quartiles = list(np.percentile( arr, [25, 50, 75] ))
            quartiles = [arr.min()] + quartiles + [arr.max()]
            whole_proteome_percentage = concatenated_lcd_lengths[aa] / total_proteome_length * 100

            output.write('\t'.join( [proteome, organism, aa] + [str(round(x, 3)) for x in quartiles] + [str( np.mean(arr) ), str(whole_proteome_percentage)] ) + '\n')
    
    output.close()
    
def calculate_perc_lcd(proteome, prot_df, lcd_content_df, seqs_df, concatenated_lcd_lengths):
    
    for aa in amino_acids:
        h = open(proteome + '.fasta_' + aa + '_RESULTS.tsv')    #NEED TO CHANGE BACK
        for i in range(7):
            h.readline()
        header = h.readline()
        
        perc_lcd_content = 0
        for line in h:
            items = line.rstrip().split('\t')
            temp = items[0].split('|')
            prot = temp[1]
            domain_bounds = items[2].replace('(', '').replace(')', '').split('-')
            lcd_len = int(domain_bounds[1]) - int(domain_bounds[0]) + 1
            perc_lcd_content += lcd_len / len(seqs_df[prot]) * 100
            prot_df[prot][aa] += lcd_len / len(seqs_df[prot]) * 100
            concatenated_lcd_lengths[aa] += lcd_len
        h.close()
            
            
    for prot in prot_df:
        for aa in amino_acids:
            lcd_content_df[aa].append( prot_df[prot][aa] )
                
    return prot_df, lcd_content_df, concatenated_lcd_lengths
                    
                    
def get_prots(proteome):

    h = open(proteome + '.fasta')
    df = {}
    prots = []
    
    for seq_record in SeqIO.parse(h, 'fasta'):
        id = str(seq_record.id)
        items = id.split('|')
        uniprot = items[1]
        seq = str(seq_record.seq)
        df[uniprot] = seq
        prots.append(uniprot)
        
    h.close()
        
    return prots, df


def get_proteomes():

    h = open('Viruses.txt')
    proteomes = []
    for line in h:
        items = line.rstrip().split(' ')
        filename = items[-1]
        if filename.endswith('.fasta.gz') and 'DNA' not in filename:
            accession, junk, extension = filename.split('.')
            proteomes.append( accession )
    h.close()
    
    return proteomes


if __name__ == '__main__':
    main()