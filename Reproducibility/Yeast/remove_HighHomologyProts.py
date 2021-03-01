
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'

def main():
    
    high_homol_df = get_high_homol_prots()
    
    for aa in amino_acids:
        for res2 in amino_acids.replace(aa, ''):
            output = open('Scerevisiae_' + aa + '-primary_' + res2 + '-secondary_LCD-containing_proteins_80percHomologyProtsRemoved.txt', 'w')
            try:
                h = open('Scerevisiae_' + aa + '-primary_' + res2 + '-secondary_LCD-containing_proteins.txt')
            except:
                continue
            
            for line in h:
                prot = line.rstrip()
                if prot not in high_homol_df[aa][res2]:
                    output.write(prot + '\n')
                    
            output.close()
            h.close()

def get_high_homol_prots():
    
    df = {aa:{res:[] for res in amino_acids} for aa in amino_acids}
    
    h = open('HighHomologyProts_80percentThreshold_SecondaryAA_Prots.tsv')
    header = h.readline()
    
    for line in h:
        aa, res2, prot = line.rstrip().split('\t')
        items = prot.split('_')
        uniprot = items[0]
        df[aa][res2].append(uniprot)
        
    h.close()
    
    return df
    



if __name__ == '__main__':
    main()