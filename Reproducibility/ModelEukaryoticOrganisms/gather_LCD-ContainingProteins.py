

def main():
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    proteomes = ['UP000002311_Scerevisiae_NoIsoforms', 'UP000001940_6239_Celegans_NoIsoforms', 'UP000000803_7227_Dmelanogaster_NoIsoforms', 'UP000000437_7955_Drerio_NoIsoforms', 'UP000186698_8355_Xlaevis_NoIsoforms', 'UP000000589_10090_Mmusculus_NoIsoforms', 'UP000005640_9606_Hsapiens_NoIsoforms']
    abbrevs = ['Scerevisiae', 'Celegans', 'Dmelanogaster', 'Drerio', 'Xlaevis', 'Mmusculus', 'Hsapiens']
    
    index = 0
    for proteome in proteomes:
        for aa in amino_acids:
            h = open(proteome + '_' + aa + '_LCD-Composer_RESULTS.tsv')
            for i in range(7):
                h.readline()
                
            output = open(abbrevs[index] + '_' + aa + '_LCD-containing_proteins.txt', 'w')
                
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