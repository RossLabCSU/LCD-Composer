
def main():

    output = open('Clustered_Benchmark_Sequences.FASTA', 'w')
    
    pep_num = 1
    for num_muts in range(1, 20):
    
        num_x = 20-num_muts
        seq = 'G'*int(num_x/2) + 'X'*num_muts
        while len(seq) < 20:
            seq += 'G'
        
        seq = 'G'*40 + seq + 'G'*40
        output.write('>PepNum' + str(pep_num) + '_NumMuts' + str(num_muts) + '\n')
        output.write(seq + '\n')
        
        pep_num += 1

    output.close()

if __name__ == '__main__':
    main()