

def main():

    output = open('RUN_LCD-Composer_ClusteredBenchmarkSequences.bat', 'w')
    for lin_disp in [str(round(x/10, 1)) for x in range(0, 11)]:
        output.write('python LCD-Composer_Modified_Output_for_ClusteredBenchmarkSequences.py Clustered_Benchmark_Sequences.FASTA Clustered_BenchmarkSequences_LCD-Composer_RESULTS_LinDisp_' + lin_disp.replace('.', 'p') + ' -a G -d ' + lin_disp + ' -i 100\n')
        
    output.close()
        
    
if __name__ == '__main__':
    main()