
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.lines import Line2D

def main():

    h = open('20aaBenchmarkSequences_LCD-Composer_RESULTS.tsv')
    for i in range(7):
        h.readline()
    header = h.readline()
    
    vals = {}
    for line in h:
        items = line.rstrip().split('\t')
        id, seq, boundaries, perc_comp, normed_stddev, *remainder = items
        num_muts = seq.count('X')
        vals[num_muts] = vals.get(num_muts, [])
        vals[num_muts].append( float(normed_stddev) )
    h.close()

    plot_linegraph(vals)
    
        
def plot_linegraph(vals):

    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
    colors = colors[::-1]
    
    index = 0
    for num_muts in sorted(list(vals.keys())[:10], reverse=True):
        yvals = []
        for threshold in [x/10 - 0.000001 for x in range(0,11)]:
            hits = [1 if x>threshold else 0 for x in vals[num_muts]]
            count = sum(hits)
            yvals.append((1-count/len(vals[num_muts])) * 100)
        plt.plot(yvals, color=colors[index])
        index += 1
        if index >= len(colors):
            index = 0
        
    plt.xticks([x for x in range(11)], labels=[str(round(x/10, 1)) for x in range(0, 11)], fontname='Arial', fontsize=16)
    plt.yticks(fontname='Arial', fontsize=16)
    plt.xlabel('Linear Dispersion Threshold', fontname='Arial', fontsize=18)
    plt.ylabel('Percentage of Sequences\nExcluded', fontname='Arial', fontsize=18)
    
    labels = []
    for x in range(1, 10):
        labels.append(str(x) + ', ' + str(20-x))

    labels.append('10')

    lines = [Line2D([0], [0], color=c, linewidth=3) for c in colors[::-1]]

    leg = plt.legend(lines, labels, prop={'family':'Arial', 'size':12.5})
    leg.set_title("Number of\nMutations", prop={'family':'Arial', 'size':14})
    
    plt.tight_layout()
    fig = plt.gcf()
    fig.set_size_inches(7, 5)
    plt.savefig('Fig S2 - Discriminatory Effect of Linear Dispersion_20aaBenchmarkSequences.tiff', bbox_inches='tight', dpi=600)
    plt.close()
    
    
if __name__ == '__main__':
    main()