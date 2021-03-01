
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'

def main():

    email = input('\nEnter your email address (required for ClustalW server): ')
    output = open('Run_Yeast_PrimaryLCDs_Alignments_ClustalW_BatchFile.bat', 'w')
    for aa in amino_acids:
        h = open('Scerevisiae_' + aa + '_LCD-containing_PROTEIN_SEQUENCES.FASTA')
        prots = []
        for line in h:
            prots.append(line.rstrip())
            
        if len(prots) > 0:
            output.write('python .\ClustalW_API.py --email ' + email + ' --stype protein --sequence Scerevisiae_' + aa + '_LCD-containing_PROTEIN_SEQUENCES.FASTA --outfile ' + aa + '_PrimaryLCD_Prots --outformat pim\n')

    output.close()


if __name__ == '__main__':
    main()