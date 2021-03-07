
organisms = ['Scerevisiae', 'Celegans', 'Dmelanogaster', 'Drerio', 'Xlaevis', 'Mmusculus', 'Hsapiens']
proteomes = ['UP000002311_Scerevisiae_NoIsoforms', 'UP000001940_6239_Celegans_NoIsoforms', 'UP000000803_7227_Dmelanogaster_NoIsoforms', 'UP000000437_7955_Drerio_NoIsoforms', 'UP000186698_8355_Xlaevis_NoIsoforms', 'UP000000589_10090_Mmusculus_NoIsoforms', 'UP000005640_9606_Hsapiens_NoIsoforms']
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'

def main():

    output = open('Run_LCD-Composer_ModelEukaryoticOrganisms_PrimaryAA-SecondaryAA.bat', 'w')

    for i in range(len(organisms)):
        organism = organisms[i]
        file = proteomes[i]
        for aa in amino_acids:
            for secondary_aa in amino_acids.replace(aa, ''):
                results_file = organism + '_' + aa + '-primary_' + secondary_aa + '-secondary_LCD-Composer_RESULTS'
                output.write('python LCD-Composer.py ' + file + '.fasta ' + results_file + ' -a ' + aa + '_' + secondary_aa + ' -c 40_20\n')
                    
    output.close()
    

if __name__ == '__main__':
    main()