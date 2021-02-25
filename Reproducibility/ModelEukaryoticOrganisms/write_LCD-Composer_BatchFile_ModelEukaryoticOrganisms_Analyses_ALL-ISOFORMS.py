
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'

def main():

    output = open('Run_LCD-Composer_on_ModelEukaryotes_AAvaried_ALL-ISOFORMS.bat', 'w')
    proteomes = ['UP000002311_Scerevisiae_IsoformsIncluded', 'UP000001940_Celegans_IsoformsIncluded', 'UP000000803_Dmelanogaster_IsoformsIncluded', 'UP000000437_Drerio_IsoformsIncluded', 'UP000186698_Xlaevis_IsoformsIncluded', 'UP000000589_Mmusculus_IsoformsIncluded', 'UP000005640_Hsapiens_IsoformsIncluded']

    for proteome in proteomes:
        for aa in amino_acids:
            results_file = proteome + '_' + aa + '_LCD-Composer_RESULTS'
            output.write('python LCD-Composer.py ' + proteome + '.fasta ' + results_file + ' -a ' + aa + '\n')

    output.close()
            
if __name__ == '__main__':
    main()