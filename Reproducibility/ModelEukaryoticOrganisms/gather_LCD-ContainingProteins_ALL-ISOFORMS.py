

def main():
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    proteomes = ['UP000002311_Scerevisiae_IsoformsIncluded', 'UP000001940_Celegans_IsoformsIncluded', 'UP000000803_Dmelanogaster_IsoformsIncluded', 'UP000000437_Drerio_IsoformsIncluded', 'UP000186698_Xlaevis_IsoformsIncluded', 'UP000000589_Mmusculus_IsoformsIncluded', 'UP000005640_Hsapiens_IsoformsIncluded']
    abbrevs = ['Scerevisiae', 'Celegans', 'Dmelanogaster', 'Drerio', 'Xlaevis', 'Mmusculus', 'Hsapiens']
    
    index = 0
    for proteome in proteomes:
        for aa in amino_acids:
            h = open(proteome + '_' + aa + '_LCD-Composer_RESULTS.tsv')
            for i in range(7):
                h.readline()
                
            output = open(abbrevs[index] + '_' + aa + '_LCD-containing_proteins_ALL-ISOFORMS.txt', 'w')
                
            prots = set()
            header = h.readline()
            for line in h:
                id, seq, boundaries, final_comp, final_stddev, *remainder = line.rstrip().split('\t')
                junk, uniprot_id, junk = id.split('|')
                prots.add( uniprot_id )
            h.close()
            
            for prot in prots:
                output.write(prot + '\n')
            output.close()
                
        index += 1

if __name__ == '__main__':
    main()