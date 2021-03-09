
from Bio import SeqIO
import pickle
import random
from goatools import obo_parser
from goatools.go_enrichment import GOEnrichmentStudy

proteomes = ['UP000002311_Scerevisiae_NoIsoforms', 'UP000001940_6239_Celegans_NoIsoforms', 'UP000000803_7227_Dmelanogaster_NoIsoforms', 'UP000000437_7955_Drerio_NoIsoforms', 'UP000186698_8355_Xlaevis_NoIsoforms', 'UP000000589_10090_Mmusculus_NoIsoforms', 'UP000005640_9606_Hsapiens_NoIsoforms']
organisms = ['Scerevisiae', 'Celegans', 'Dmelanogaster', 'Drerio', 'Xlaevis', 'Mmusculus', 'Hsapiens']
pop_sizes = [6049, 19818, 13806, 25698, 43236, 21989, 20600]    # SIZES OF RESPECTIVE PROTEOMES
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'

goterm_to_assoc_prots = pickle.load(open('GOterms_dict_with_ActualAssociatedProtSet_as_Keys.dat', 'rb'))


def main():

    sampsize_df = get_LCD_sampsizes()
    full_proteome_lists = get_proteomes()
    output = open('Results_from_PvalThreshold_Testing_for_GOtermSampling.tsv', 'w')
    output.write('\t'.join( ['# GO', 'NS' 'enrichment', 'name', 'ratio_in_study', 'ratio_in_pop', 'p_uncorrected', 'depth', 'study_count', 'p_bonferroni', 'p_sidak', 'p_holm', 'study_items'] ) + '\n')
    go_obo = 'go-basic.obo'
    go = obo_parser.GODag(go_obo)
    
    pval_thresh_df = {}
    for organism in organisms:
        goterm_to_assoc_prots_sampling_list = list(goterm_to_assoc_prots[organism].items())   #create list of dictionary items for random.choice
        pval_thresh_df[organism] = {}
        background_prots = full_proteome_lists[organism]
        assoc = get_assoc(organism)
        for aa in amino_acids:
            pval_thresh_df, output = find_min_pval_thresh(organism, aa, background_prots, sampsize_df, pval_thresh_df, output, goterm_to_assoc_prots_sampling_list, assoc, go)
            
    output.close()
    
    pf = open('OrganismSpecific_and_LCDspecific_NumAssocProteins_Thresholds_for_GOsignificance.dat', 'wb')
    pickle.dump(pval_thresh_df, pf)
    pf.close()
  
    
def get_LCD_sampsizes():
    
    df = {}
    for organism in organisms:
        df[organism] = {}
        for aa in amino_acids:
            samp_size = 0
            h = open(organism + '_' + aa + '_LCD-containing_proteins.txt')
            for line in h:
                prot = line.rstrip()
                if prot != '':
                    samp_size += 1
            h.close()
            
            df[organism][aa] = samp_size

    return df
    
    
def find_min_pval_thresh(organism, aa, background_prots, sampsize_df, pval_thresh_df, output, goterm_to_assoc_prots_sampling_list, assoc, go):
    
    num_in_study = 1
    sidak_pval = 1.0
    lcd_prots = get_hits(organism, aa)
    
    # EXITS FUNCTION AND SKIPS CALCULATIONS FOR GROUPS WITH NO PROTEINS
    if len(lcd_prots) == 0:
        return pval_thresh_df, output
        
    while sidak_pval > 0.05 - 0.0000001:
        target_go_id, go_prot_set = get_sampled_GOterm( organism, num_in_study, goterm_to_assoc_prots_sampling_list )
        modified_lcd_prots = modify_lcd_protlist( lcd_prots, go_prot_set )

        methods = ["bonferroni", "sidak", "holm", "fdr"]
        g = GOEnrichmentStudy(background_prots, assoc, go, propagate_counts=False, alpha=0.05)
        g_res = g.run_study(modified_lcd_prots)
        for rec in g_res:
            rec = str(rec)
            items = rec.split('\t')
            go_id = items[0]

            if go_id != target_go_id:
                continue

            enr_or_pur = items[2]
            p_uncorr = float(items[6])
            sidak_pval = float(items[10])
                
            output.write(organism + '\t' + aa + '\t' + rec + '\n')
            
        num_in_study += 1
            
    pval_thresh_df[organism][aa] = num_in_study - 1     # THE -1 HERE ACCOUNTS FOR THE FACT THAT +1 WAS ADDED TO num_in_study FOR THE LAST WHILE LOOP ITERATION EVEN THOUGH YOU WANT TO STORE THE VALUE CORRESPONDING TO THAT ITERATION.
            
    return pval_thresh_df, output
                
                
def modify_lcd_protlist( lcd_prots, go_prot_set ):
    """Function that modifies the lcd_prot list by replacing randomly selected proteins with 
    the proteins from the GO term set. Effectively plants the proteins of interest within the LCD set
    so that the all of the GO term associated proteins will be identified (simulates maximum enrichment
    for that set of proteins).
    """
    index = 0
    for prot in go_prot_set:
        if prot in lcd_prots:
            continue
            
        replaced_prot = random.choice(lcd_prots)
        while replaced_prot in go_prot_set:
            replaced_prot = random.choice(lcd_prots)
            
        index = lcd_prots.index(replaced_prot)
        lcd_prots = lcd_prots[:index] + [prot] + lcd_prots[index+1:]
            
    return lcd_prots
           
                
def get_sampled_GOterm(organism, num_in_study, goterm_to_assoc_prots_sampling_list):

    # LOOP SHOULD RUN UNTIL A GO TERM WITH THE PROPER NUMBER OF PROTEINS IS SAMPLED. THEN LOOP SHOULD EXIT AND THE FINAL go_id AND prot_set VALUES SHOULD CORRESPOND TO A VALID SAMPLED GO TERM WITH THE PROPER NUMBER OF ASSOCIATED PROTEINS.
    num_prots = 0
    while num_prots != num_in_study:
        go_id, prot_set = random.choice(goterm_to_assoc_prots_sampling_list)
        num_prots = len(prot_set)
        
    return go_id, prot_set
                
    
def get_proteomes():

    all_proteins = {}
    for i in range(len(organisms)):
        h = open(proteomes[i] + '.fasta')
        prots = []
        for seq_record in SeqIO.parse(h, 'fasta'):
            id = str(seq_record.id)
            seq = str(seq_record.seq)
            if seq[-1] == '*':
                seq = seq[:-1]
            junk, uniprot, junk = id.split('|')
            prots.append( uniprot )

        organism = organisms[i]
        all_proteins[organism] = prots

    return all_proteins
    
    
def get_hits(proteome, aa):

    prots = []
    h = open(proteome + '_' + aa + '_LCD-containing_proteins.txt')
    for line in h:
        prots.append(line.rstrip())
        
    return prots
    
    
def get_assoc(abbrev):
    """
    Reads the s_cerevisiae gene association file (GAF).
    Return is a dictionary with SGD ORFs as keys, and a set of associated GO terms as values (this is a formal Python set, not a list)
    
    Each list represents a separate GO term associated with the gene key, and contains the following:
    0) Database ("SGD")
    1) SGD_ID
    2) Common gene name (e.g. SUP35)
    3) Qualifier (e.g. 'NOT', 'contributes_to', 'colocalizes_with', etc. This is optional)
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
    14) Assigned by (always 'SGD' for this file)
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
