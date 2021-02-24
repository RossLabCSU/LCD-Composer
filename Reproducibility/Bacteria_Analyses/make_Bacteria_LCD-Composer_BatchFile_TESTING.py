
import os
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'

def main():

    output = open('RUN_LCD-Composer_Bacteria_Batch.bat', 'w')

    for (dirname, dir, files) in os.walk('.'):
        for file in files:
            if file.endswith('.fasta'):
                uniprot_name, extension = file.split('.')
                for aa in amino_acids:
                    output.write('python LCD-Composer.py ' + file + ' ' + uniprot_name + '_' + aa + '_RESULTS -a ' + aa + '\n')
    output.close()

if __name__ == '__main__':
    main()