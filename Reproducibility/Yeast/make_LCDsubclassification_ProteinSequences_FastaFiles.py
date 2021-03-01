
from Bio import SeqIO
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'

def main():

    seqs = get_proteome()
    for aa in amino_acids:
        for sec_aa in amino_acids.replace(aa, ''):
            try:
                h = open('Scerevisiae_' + aa + '-primary_' + sec_aa + '-secondary_LCD-containing_proteins.txt')
            except:
                continue
            output = open('Scerevisiae_' + aa + '-primary_' + sec_aa + '-secondary_PROTEIN_SEQUENCES.FASTA', 'w')
            for line in h:
                id = line.rstrip()
                common_name, seq = seqs[id]
                output.write('>' + id + '_' + common_name + '\n')
                output.write(seq + '\n')
            h.close() 
            output.close()

def get_proteome():

    h = open('UP000002311_Scerevisiae_NoIsoforms.FASTA')
    df = {}
    full_ids = {}
    for seq_record in SeqIO.parse(h, 'fasta'):
        id = str(seq_record.id)
        desc = str(seq_record.description)
        items = desc.split(' ')
        for item in items:
            if 'GN=' in item:
                junk, gene_name = item.split('=')
        junk, uniprot, junk = id.split('|')
        seq = str(seq_record.seq)
        
        if seq[-1] == '*':
            seq = seq[:-1]

        name = uniprot + '_' + gene_name
        df[uniprot] = (gene_name, seq)
        
    return df


if __name__ == '__main__':
    main()
