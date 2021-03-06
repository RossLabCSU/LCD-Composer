
proteomes = ['UP000002311_Scerevisiae_NoIsoforms', 'UP000001940_6239_Celegans_NoIsoforms', 'UP000000803_7227_Dmelanogaster_NoIsoforms', 'UP000000437_7955_Drerio_NoIsoforms', 'UP000186698_8355_Xlaevis_NoIsoforms', 'UP000000589_10090_Mmusculus_NoIsoforms', 'UP000005640_9606_Hsapiens_NoIsoforms']
organisms = ['Scerevisiae', 'Celegans', 'Dmelanogaster', 'Drerio', 'Xlaevis', 'Mmusculus', 'Hsapiens']
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'

def main():

    output = open('TableS3_ModelEukaryoticOrganisms_All_LCDs.tsv', 'w')
    output.write('\t'.join( ['Organism', 'LCD Class', 'Protein ID', 'Domain Boundaries', 'AA Used in Search', 'Final Comp of Search AA', 'Final Dispersion of Search AA'] ) + '\n')

    for i in range(len(organisms)):
        proteome = proteomes[i]
        organism = organisms[i]
        organism = organism[0] + '. ' + organism[1:]
        for aa in amino_acids:
            h = open(proteome + '_' + aa + '_LCD-Composer_RESULTS.tsv')
            for j in range(7):
                h.readline()
            header = h.readline()
            
            for line in h:
                items = line.rstrip().split('\t')
                temp, uniprot, *remainder = items[0].split('|')
                output.write('\t'.join( [organism, aa, uniprot] + items[1:3] + [aa] + items[3:5] ) + '\n')
                
            h.close()
            
    output.close()

if __name__ == '__main__':
    main()