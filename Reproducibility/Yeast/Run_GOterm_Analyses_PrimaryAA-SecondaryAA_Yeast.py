
from goatools import obo_parser
from goatools.go_enrichment import GOEnrichmentStudy
from Bio import SeqIO
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'

def main():

    go_obo = 'go-basic.obo'
    go = obo_parser.GODag(go_obo)
    
    organism = 'Scerevisiae'
    proteome = 'UP000002311_Scerevisiae_NoIsoforms'
    
    background_prots = get_background_prots(proteome)
    assoc = get_assoc(organism)
    for aa in amino_acids:
        for secondary_aa in amino_acids.replace(aa, ''):
        
            lcd_prots = get_hits(organism, aa, secondary_aa)
            
            if len(lcd_prots) > 0:
                methods = ["bonferroni", "sidak", "holm", "fdr"]
                g = GOEnrichmentStudy(background_prots, assoc, go, propagate_counts=False, alpha=0.05)
                g_res = g.run_study(lcd_prots)
                file_name = organism + '_' + aa + '-primary_' + secondary_aa + '-secondary_GO_RESULTS.tsv'
                g.wr_tsv(file_name, g_res)
           
    
def get_hits(proteome, aa, secondary_aa):

    prots = []
    h = open(proteome + '_' + aa + '-primary_' + secondary_aa + '-secondary_LCD-containing_proteins.txt')
    for line in h:
        prots.append(line.rstrip())
    h.close()
        
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