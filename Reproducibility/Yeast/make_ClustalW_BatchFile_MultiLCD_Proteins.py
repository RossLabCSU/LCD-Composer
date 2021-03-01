
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'

def main():

    email = input('\nEnter your email address (required for ClustalW server): ')
    output = open('Run_Yeast_MultiLCD_Alignments_ClustalW_BatchFile.bat', 'w')
    for aa in amino_acids:
        for sec_aa in amino_acids.replace(aa, ''):
            try:
                h = open('Scerevisiae_' + aa + '_' + sec_aa + '_MultiLCD_PROTEIN_SEQUENCES.FASTA')
            except:
                continue
            prots = []
            for line in h:
                prots.append(line.rstrip())
                
            if len(prots) > 0:
                output.write('python .\ClustalW_API.py --email ' + email + ' --stype protein --sequence Scerevisiae_' + aa + '_' + sec_aa + '_MultiLCD_PROTEIN_SEQUENCES.FASTA --outfile ' + aa + '_' + sec_aa + '_MultiLCD_Prots --outformat pim\n')

    output.close()


if __name__ == '__main__':
    main()