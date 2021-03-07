
organisms = ['Scerevisiae', 'Celegans', 'Dmelanogaster', 'Drerio', 'Xlaevis', 'Mmusculus', 'Hsapiens']

def main():

    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    for organism in organisms:
        for aa in amino_acids:
            for sec_aa in amino_acids.replace(aa, ''):
                h = open(organism + '_' + aa + '-primary_' + sec_aa + '-secondary_LCD-Composer_RESULTS.tsv')
                output = open(organism + '_' + aa + '-primary_' + sec_aa + '-secondary_LCD-containing_proteins.txt', 'w')
                
                for i in range(7):
                    h.readline()
                header = h.readline()
                    
                prots = set()
                for line in h:
                    id, seq, boundaries, final_comp, final_lindisp, *remainder = line.rstrip().split('\t')
                    junk, uniprot_id, junk = id.split('|')
                    prots.add( uniprot_id )
                h.close()
                
                for prot in prots:
                    output.write(prot + '\n')
                    
                output.close()


if __name__ == '__main__':
    main()