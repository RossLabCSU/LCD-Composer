

def main(args):

    seq_len = args.sequence_length
    mut_range = args.mutation_range
    
    # CHECK FOR VALID mutation_range
    if mut_range < 1 or mut_range > seq_len-1:
        print('Invalid mutation_range. This must be an integer between 1 and the specified sequence length minus 1 (inclusive).')
        exit()
        
    proceed = input('\nWith default parameters, this will generate all possible 20aa sequence combinations which results in a file size of ~0.15GB.\n\nFile sizes become quickly extremely large as sequence length is increased. For example, the file generated using a sequence length =30 and mutation range =15 will contain all necessary 30aa sequences to fully test the linear dispersion parameter, but the file is >35GB in size.\n\nThe mutation range =15 only generates sequences up to 50% composition for the two amino acids, which is the minimum necessary to test the linear dispersion parameter since reciprocal sequences will result in identical linear dispersion values. However, if you wish to generate the complete set of possible 30aa sequences use the "-r 29" flag, and the resulting file size should be >62GB.\n\nDo you wish to proceed? Enter yes or no:  ')

    # CHECK FOR VALID RESPONSE
    if proceed.lower() not in ['yes', 'y', 'no', 'n']:
        print('Invalid response. You must respond by answering "yes" or "no".')
        exit()
        
    if proceed.lower() == 'yes' or proceed.lower() == 'y':
            
        from itertools import combinations
        output = open(str(seq_len) + 'aaBenchmarkSequences.FASTA', 'w')

        # seq_len = 20
        pep_num = 0
        for num_muts in range(1, mut_range+1):
            print("Generating sequences with " + str(num_muts) + " number of mutations (X's)", str(datetime.datetime.now()))
            for position_comb in combinations([x for x in range(seq_len)], num_muts):
                pep_num += 1
                seq = 'G'*seq_len
                for pos in position_comb:
                    seq = seq[:pos] + 'X' + seq[pos+1:]
                
                output.write('>PepNum' + str(pep_num) + '_NumMuts' + str(num_muts) + '\n')
                output.write(seq + '\n')
                
        output.close()
        
    else:
        print('Exiting program without running...')
        exit()
        
        
    
def get_args(arguments):
    parser = argparse.ArgumentParser(description='Generate all possible Benchmark sequences consisting of only 2 types of amino acid symbols', prog='generateBenchmarkSequences')

    parser.add_argument('-r', '--mutation_range', type=int, default=19, 
                        help="""Mutation range (integer). Represents the maximum number of "mutations" that you wish to generate sequences for.
                                    Benchmark sequences are 30aa sequences starting with 29 G residues and 1 X residue (representing 1 "mutation").
                                    All possible arrangements of G+X are generated with these compositions, then the number of "mutations" is 
                                    increased up to (and including) the number specified by --mutation_range.
                                    
                                    It is recommended that you only do a maximum of 15, since the linear dispersion for all subsequent sequences
                                    will be identical to the reciprocal sequence which will already be analyzed with the mutations 1-15.
                                    e.g. the linear dispersion for the sequence GGGXGXGGG is identical to the linear dispersion for the sequence XXXGXGXXX.""")
                                    
    parser.add_argument('-l', '--sequence_length', type=int, default=20, 
                        help="""Mutation range (integer). Represents the maximum number of "mutations" that you wish to generate sequences for.
                                    e.g. Benchmark sequences that are 20aa in length will start with sequences consisting of 19 G residues and 1 X residue (representing 1 "mutation").
                                    All possible arrangements of G+X are generated with these compositions, then the number of "mutations" is 
                                    increased up to (and including) the number specified by --mutation_range.
                                    
                                    For sequence lengths larger than 20aa, it is recommended that you only do a mutation range up to 1/2 the sequence length, 
                                    since the linear dispersion for all subsequent sequences will be identical to the reciprocal sequence, which will 
                                    already be analyzed for a smaller number of mutations.
                                    e.g. the linear dispersion for the sequence GGGXGXGGG is identical to the linear dispersion for the sequence XXXGXGXXX.""")
                                    
    args = parser.parse_args(arguments)
    
    return args

            
if __name__ == '__main__':
    import sys, argparse, datetime
    args = get_args(sys.argv[1:])
    main(args)
