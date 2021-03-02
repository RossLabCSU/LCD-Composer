
from Bio import SeqIO
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
import datetime
import pickle

def main():

    df = {}
    win_size = 20
    print('Start time:', str(datetime.datetime.now()))
    for prim_aa in amino_acids:
        print(prim_aa, 'Start time:', str(datetime.datetime.now()))
        secondary_aas = amino_acids.replace(prim_aa, '')
        df[prim_aa] = {}
        for secondary_aa in secondary_aas:
            df[prim_aa][secondary_aa] = 0
        total_windows = 0

        h = open('UP000002311_Scerevisiae_NoIsoforms.fasta')
        for seq_record in SeqIO.parse(h, 'fasta'):
            seq = str(seq_record.seq)
            if seq[-1] == '*':
                seq = seq[:-1]
                
            for i in range(len(seq) - win_size):
                window = seq[i:i+win_size]
                counts = [window.count(aa) for aa in secondary_aas]
                max_count = max(counts)
                if counts.count(max_count) > 1:
                    continue
                ind = counts.index(max_count)
                secondary_aa = secondary_aas[ind]
                df[prim_aa][secondary_aa] += 1
                total_windows += 1
        h.close()
        df[prim_aa]['Total Windows'] = total_windows

    pf = open('Background_X-rich_Window_SecondaryAA_Frequencies.dat', 'wb')
    pickle.dump(df, pf)
    pf.close()


if __name__ == '__main__':
    main()