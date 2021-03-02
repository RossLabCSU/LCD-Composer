
import matplotlib.pyplot as plt
import seaborn as sns

def main():

    xvals = []
    yvals = []
    i = 0
    labels = []
    seqs = []
    for lin_disp in [str(round(x/10, 1)) for x in range(0, 11)]:
        h = open('Clustered_BenchmarkSequences_LCD-Composer_RESULTS_LinDisp_' + lin_disp.replace('.', 'p') + '.tsv')
        
        for j in range(7):
            h.readline()
            
        header = h.readline()
        toggle = 0
        for line in h:
            items = line.rstrip().split(',')
            id, full_seq, domain_seqs, domain_bounds = line.rstrip().split('\t')
            domain_seqs = domain_seqs.split('_')
            seqs.append(full_seq[20:-20])

            if len(domain_seqs) > 1 and toggle==0:
                toggle = 1
                xvals.append(i)
                yvals.append(full_seq.count('X'))
                labels.append(lin_disp)
                
        #Catches 1.0 threshold case, which doesn't exclude any sequences and doesn't split any of the sequences
        if toggle==0:
            xvals.append(i)
            yvals.append(full_seq.count('X')+1)
            labels.append(lin_disp)
                
        h.close()
        i += 1
        
    plt.plot(xvals, yvals)
    plt.scatter(xvals, yvals)
    plt.xticks([x for x in range(0, 11)], labels=[str(round(x/10, 1)) for x in range(0, 11)], fontname='Arial', fontsize=16)
    plt.yticks([x for x in range(1, 20)], labels=[seqs[i] + ', ' + str(i+1) for i in range(0, 20)], fontname='Arial', fontsize=12)
    plt.ylim(1, 20)
    plt.xlabel('Linear Dispersion Threshold', fontname='Arial', fontsize=18)
    plt.ylabel('# of MutationsRequired\nto Split Domains', fontname='Arial', fontsize=18)
    plt.grid(color='0.7')
    fig = plt.gcf()
    fig.set_size_inches(15, 5)
    plt.tight_layout()
    plt.savefig('Fig S3 - Breakpoints_for_ClusteredBenchmarkSequences.tiff', bbox_inches='tight', dpi=600)
    plt.close()
                
            
if __name__ == '__main__':
    main()