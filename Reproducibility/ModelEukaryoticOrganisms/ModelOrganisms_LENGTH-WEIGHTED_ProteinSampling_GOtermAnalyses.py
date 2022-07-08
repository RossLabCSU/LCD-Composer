
import numpy as np
from Bio import SeqIO
from goatools import obo_parser
from goatools.go_enrichment import GOEnrichmentStudy
import datetime

amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
aa_names = {'A':'Alanine (A)', 'C':'Cysteine (C)', 'D':'Aspartic Acid (D)', 'E':'Glutamic Acid (E)',
        'F':'Phenylalanine (F)', 'G':'Glycine (G)', 'H':'Histidine (H)', 'I':'Isoleucine (I)', 
        'K':'Lysine (K)', 'L':'Leucine (L)', 'M':'Methionine (M)', 'N':'Asparagine (N)',
        'P':'Proline (P)', 'Q':'Glutamine (Q)', 'R':'Arginine (R)', 'S':'Serine (S)',
        'T':'Threonine (T)', 'V':'Valine (V)', 'W':'Tryptophan (W)','Y':'Tyrosine (Y)'}
        
proteomes = ['UP000002311_Scerevisiae_NoIsoforms', 'UP000001940_6239_Celegans_NoIsoforms', 'UP000000803_7227_Dmelanogaster_NoIsoforms', 'UP000000437_7955_Drerio_NoIsoforms', 'UP000186698_8355_Xlaevis_NoIsoforms', 'UP000000589_10090_Mmusculus_NoIsoforms', 'UP000005640_9606_Hsapiens_NoIsoforms']
abbrevs = ['Scerevisiae', 'Celegans', 'Dmelanogaster', 'Drerio', 'Xlaevis', 'Mmusculus', 'Hsapiens']

def main():

    go_obo = 'go-basic.obo'
    go = obo_parser.GODag(go_obo)

    full_proteome_lists, lw_prot_lists = get_proteomes()

    for i in range(len(abbrevs)):
        lw_prot_list = lw_prot_lists[abbrevs[i]]
        prots_output = open(abbrevs[i] + '_RandomlySelectedProts_All_Iterations_LENGTH_WEIGHTED.csv', 'w')
        prots_output.write(','.join( ['LCD Class', 'Iteration Number', 'Sampled Proteins (underscore-delimited)'] ) + '\n')
        output = open(abbrevs[i] + '_RandomlySampledProteins_All_Iterations_GO_RESULTS_LENGTH_WEIGHTED.tsv', 'w')
        output.write('\t'.join( ['# GO', 'NS', 'enrichment', 'name', 'ratio_in_study', 'ratio_in_pop', 'p_uncorrected', 'depth', 'study_count', 'p_bonferroni', 'p_sidak', 'p_holm', 'study_items', 'Iteration Number', 'LCD Class'] ) + '\n')

        assoc = get_assoc(abbrevs[i])

        for iteration in range(1, 1001):
            print('\n\nIteration: ', iteration, str(datetime.datetime.now()))
            #RUN GO term analysis
            for aa in amino_acids:
            
                lcd_prots = get_hits(abbrevs[i], aa)
                num_prots = len(lcd_prots)
                background_prots = full_proteome_lists[ abbrevs[i] ]
                
                random_prots = length_weighted_prot_sample(lw_prot_list, num_prots)

                if len(random_prots) > 0:
                    methods = ["bonferroni", "sidak", "holm", "fdr"]
                    g = GOEnrichmentStudy(background_prots, assoc, go, propagate_counts=False, alpha=0.05)
                    g_res = g.run_study(random_prots)
                    for rec in g_res:
                        rec = str(rec)
                        items = rec.split('\t')
                        enr_or_pur = items[2]
                        p_uncorr = float(items[6])
                        p_sidak = float(items[10])
                        if enr_or_pur == 'e' and p_sidak < 0.05+0.000000000001:
                            output.write(rec + '\t' + str(iteration) + '\t' + aa + '\n')
                    
                prots_output.write(aa + ',' + str(iteration) + ',' + '_'.join(random_prots) + '\n')
        prots_output.close()
            

def length_weighted_prot_sample(lw_prot_list, num_prots):
    
    sampled_prots = []
    while len(sampled_prots) < num_prots:
        prot_ind = np.random.choice(len(lw_prot_list))
        prot = lw_prot_list[prot_ind]
        if prot not in sampled_prots:
            sampled_prots.append( prot )
            
    return sampled_prots


def get_proteomes():

    all_proteins = {}
    lw_prot_lists = {}
    for i in range(len(abbrevs)):
        h = open(proteomes[i] + '.fasta')
        prots = []
        lw_prot_list = []
        for seq_record in SeqIO.parse(h, 'fasta'):
            id = str(seq_record.id)
            seq = str(seq_record.seq)
            if seq[-1] == '*':
                seq = seq[:-1]
            junk, uniprot, junk = id.split('|')
            prots.append( uniprot )
            for j in range(len(seq)):
                lw_prot_list.append(uniprot)

        organism = abbrevs[i]
        all_proteins[organism] = prots
        lw_prot_lists[organism] = lw_prot_list
        
    return all_proteins, lw_prot_lists
    
    
def get_hits(proteome, aa):

    prots = []
    h = open(proteome + '_' + aa + '_LCD-containing_proteins.txt')
    for line in h:
        prots.append(line.rstrip())
        
    return prots
    
        
def get_background_prots(proteome):
    
    prots = []
    h = open(proteome + '.fasta')
    for seq_record in SeqIO.parse(h, 'fasta'):
        id = str(seq_record.id)
        junk, uniprot_id, junk = id.split('|')
        prots.append( uniprot_id )
        
    return prots
    
    
def get_assoc(abbrev):
    """
    Reads the proper gene association file (GAF).
    Return is a dictionary with Protein IDs as keys, and a set of associated GO terms as values (this is a formal Python set, not a list)
    
    Each list represents a separate GO term associated with the gene key, and contains the following:
    0) Database ("UniProtKB")
    1) Protein ID
    2) Common gene name
    3) Qualifier (e.g. 'NOT', 'contributes_to', 'colocalizes_with', etc. This is optional and is sometimes an empty string)
    4) GO_ID (i.e. the GO term associated with the gene. NOTE: multiple GO terms associated with a single gene will be in separate entries in the dictionary)
    5) Literature reference that the annotation was derived from
    6) Evidence code
    7) Mystery - don't really know what this column represents and doesn't seem to match README file description on GO website
    8) Aspect (i.e. molecular function, biological process, or cellular component...indicated with abbreviation)
    9) Object symbol (e.g. 'Mitochondrial 21s RNA)
    10) Object synonym
    11) Object type (i.e. gene, protein, etc.)
    12) Taxon ID (559292 for s cerevisiae)
    13) Date of annotation
    14) Assigned by (a database)
    """

    h = open(abbrev + '.gaf')
    
    all_gos = []
    assoc = {}
    
    #creates 'pop' list and 'assoc' dictionary to pass to GOEnrichmentStudy() module in goatools
    for line in h:
        if line.startswith('!'):
            continue
        items = line.rstrip().split('\t')
        uniprot_id = items[1]
        assoc[uniprot_id] = assoc.get(uniprot_id, [])
        assoc[uniprot_id].append(items[4])
        all_gos.append(items[4])
        
    set_gos = set(all_gos)
    
    #convert associated GO term list to a set to consolidate redundant GO terms associated with a single gene
    for key in assoc:
        assoc[key] = set(assoc[key])
    
    h.close()

    return assoc

if __name__ == '__main__':
    main()
