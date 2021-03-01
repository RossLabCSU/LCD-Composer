
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'

def main():

    email = input('\nEnter your email address (required for ClustalW server): ')
    output = open('Run_Yeast_LCDsubclassification_Alignments_ClustalW_BatchFile.bat', 'w')
    for aa in amino_acids:
        for sec_aa in amino_acids.replace(aa, ''):
            h = open('Scerevisiae_' + aa + '-primary_' + sec_aa + '-secondary_PROTEIN_SEQUENCES.FASTA')
            prots = []
            for line in h:
                prots.append(line.rstrip())
                
            if len(prots) > 0:
                output.write('python .\ClustalW_API.py --email ' + email + ' --stype protein --sequence Scerevisiae_' + aa + '-primary_' + sec_aa + '-secondary_PROTEIN_SEQUENCES.FASTA --outfile ' + aa + '-primary_' + sec_aa + '-secondary_Prots --outformat pim\n')

    output.close()


if __name__ == '__main__':
    main()