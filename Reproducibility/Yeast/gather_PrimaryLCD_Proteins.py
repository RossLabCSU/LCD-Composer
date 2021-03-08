

def main():

    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    for aa in amino_acids:
        h = open('Scerevisiae_' + aa + '_40-CompositionThreshold_20-WindowSize_0p5-LinearDispersionThreshold_LCD-Composer_RESULTS.tsv')
        output = open('Scerevisiae_' + aa + '_LCD-containing_proteins.txt', 'w')
        
        for i in range(7):
            h.readline()
        header = h.readline()
            
        prots = set()
        for line in h:
            id, seq, boundaries, final_comp, final_lindisp, *remainder = line.rstrip().split('\t')
            junk, uniprot_id, junk = id.split('|')
            prots.add( uniprot_id )
        h.close()
        
        for prot in prots:
            output.write(prot + '\n')
            
        output.close()


if __name__ == '__main__':
    main()