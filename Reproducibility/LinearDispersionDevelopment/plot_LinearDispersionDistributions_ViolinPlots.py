
import matplotlib.pyplot as plt
import seaborn as sns

def main():

    h = open('20aaBenchmarkSequences_LCD-Composer_RESULTS.tsv')
    for i in range(7):
        h.readline()
    header = h.readline()

    df = {'Number of Mutations':[],
          'Linear Dispersion':[]}

    for line in h:
        items = line.rstrip().split('\t')
        id, seq, domain_bounds, final_comp, linear_disp, *remainder = items
        
        num_muts = seq.count('X')
        
        df['Number of Mutations'].append( num_muts )
        df['Linear Dispersion'].append( float(linear_disp) )

    h.close()
    
    boxplots(df)
    
    
def boxplots(df):

    sns.violinplot(x=df['Number of Mutations'], y=df['Linear Dispersion'], showfliers=False, color='#1f77b4', cut=0, bw=0.2)
    
    plt.xticks(fontname='Arial', fontsize=18)
    plt.yticks(fontname='Arial', fontsize=18)
    plt.xlabel('Number of Mutations', fontname='Arial', fontsize=20)
    plt.ylabel('Linear Dispersion', fontname='Arial', fontsize=20)
    fig = plt.gcf()
    fig.set_size_inches(15, 5)
    plt.savefig('Fig S1 - Linear Dispersion Distributions_20aaBenchmarkSequences.tiff', bbox_inches='tight', dpi=600)
    plt.close()


if __name__ == '__main__':
    main()