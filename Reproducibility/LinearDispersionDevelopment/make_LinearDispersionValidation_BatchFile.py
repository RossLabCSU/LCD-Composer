

def main():

    output = open('RUN_LinearDispersionValidation.bat', 'w')
    for seq_len in range(5, 10):
        if seq_len % 2 == 0:
            mut_range = int(seq_len/2)
        else:
            mut_range = int(seq_len/2) + 1
        output.write('python LCD-Composer_Analyze_BenchmarkSequences.py LinearDispersion_' + str(seq_len) + 'aa_Validation_RESULTS -a G -c 0 -d 0 -w ' + str(seq_len) + ' -l ' + str(seq_len) + ' -r ' + str(mut_range) + ' -m\n')
    output.close()

if __name__ == '__main__':
    main()