
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'

def main():

    output = open('Run_LCD-Composer_Yeast_PrimaryAA-SecondaryAA.bat', 'w')

    file = 'UP000002311_Scerevisiae_NoIsoforms'
    for aa in amino_acids:
        for secondary_aa in amino_acids.replace(aa, ''):
            results_file = 'Scerevisiae_' + aa + '-primary_' + secondary_aa + '-secondary_LCD-Composer_RESULTS'
            output.write('python LCD-Composer.py ' + file + '.fasta ' + results_file + ' -a ' + aa + '_' + secondary_aa + ' -c 40_20\n')
                
    output.close()
    

if __name__ == '__main__':
    main()
