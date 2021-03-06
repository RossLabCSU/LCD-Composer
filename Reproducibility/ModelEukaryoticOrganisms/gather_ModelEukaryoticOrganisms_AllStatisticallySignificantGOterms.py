
proteomes = ['UP000002311_Scerevisiae_NoIsoforms', 'UP000001940_6239_Celegans_NoIsoforms', 'UP000000803_7227_Dmelanogaster_NoIsoforms', 'UP000000437_7955_Drerio_NoIsoforms', 'UP000186698_8355_Xlaevis_NoIsoforms', 'UP000000589_10090_Mmusculus_NoIsoforms', 'UP000005640_9606_Hsapiens_NoIsoforms']
organisms = ['Scerevisiae', 'Celegans', 'Dmelanogaster', 'Drerio', 'Xlaevis', 'Mmusculus', 'Hsapiens']
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'

def main():

    output = open('TableS5_ModelEukaryoticOrganisms_AllStatisticallySignificantGOterms.tsv', 'w')

    for i in range(len(organisms)):
        proteome = proteomes[i]
        organism = organisms[i]
        for aa in amino_acids:
            try:
                h = open(organism + '_' + aa + '_GO_RESULTS.tsv')
            except:
                continue
            header = h.readline()
            if organism == 'Scerevisiae' and aa == 'A':
                header = header.rstrip().split('\t')
                output.write('\t'.join( ['Organism', 'LCD Class'] + header ) + '\n')
            
            for line in h:
                items = line.rstrip().split('\t')
                p_sidak = float(items[10])
                if p_sidak < 0.05+0.0000001:
                    name = organism[0] + '. ' + organism[1:]
                    output.write('\t'.join( [name, aa] + items ) + '\n')
                
            h.close()
            
    output.close()

if __name__ == '__main__':
    main()