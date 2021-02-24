
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'

def main():

    output = open('Run_LCD-Composer_on_ModelEukaryotes_AAvaried.bat', 'w')
    proteomes = ['UP000002311_Scerevisiae_NoIsoforms', 'UP000001940_6239_Celegans_NoIsoforms', 'UP000000803_7227_Dmelanogaster_NoIsoforms', 'UP000000437_7955_Drerio_NoIsoforms', 'UP000186698_8355_Xlaevis_NoIsoforms', 'UP000000589_10090_Mmusculus_NoIsoforms', 'UP000005640_9606_Hsapiens_NoIsoforms']

    for proteome in proteomes:
        for aa in amino_acids:
            results_file = proteome + '_' + aa + '_LCD-Composer_RESULTS'
            output.write('python LCD-Composer.py ' + proteome + '.fasta ' + results_file + ' -a ' + aa + '\n')

    output.close()
            
if __name__ == '__main__':
    main()