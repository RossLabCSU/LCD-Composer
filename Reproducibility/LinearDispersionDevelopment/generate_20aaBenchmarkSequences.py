
from itertools import combinations

def main():

    output = open('20aaBenchmarkSequences.FASTA', 'w')

    seq_len = 20
    
    n_vals = []
    pep_num = 0

    for num_muts in range(1,seq_len):
        position_combs = combinations([x for x in range(seq_len)], num_muts)
        position_combs = list(position_combs)
        n_vals.append(len(list(position_combs)))
        for position_comb in list(position_combs):
            pep_num += 1
            seq = 'G'*seq_len
            for pos in position_comb:
                seq = seq[:pos] + 'X' + seq[pos+1:]
            
            output.write('>PepNum' + str(pep_num) + '_NumMuts' + str(num_muts) + '\n')
            output.write(seq + '\n')
            
    output.close()

            
if __name__ == '__main__':
    main()
