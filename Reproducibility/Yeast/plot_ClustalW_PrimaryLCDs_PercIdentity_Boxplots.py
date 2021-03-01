
import matplotlib.pyplot as plt
import seaborn as sns
amino_acids = 'ACDEFGHIKLMNPQRSTVY'

def main():
    
    df = {'LCD Class':[],
                'Percent Identity':[]}
    for aa in amino_acids:
        aa_list = []
        h = open(aa + '_PrimaryLCD_Prots.pim.pim')
        for i in range(6):
            h.readline()
            
        for line in h:
            items = line.rstrip().split(' ')
            items = [x for x in items if x != '']
            row_num = int(items[0].replace(':', ''))
            vals = [ float(x) for x in items[2:] ]
            for i in range(0, row_num-1):
                df['LCD Class'].append(aa)
                df['Percent Identity'].append( vals[i] )
        h.close()

    boxplots(df)
            
def boxplots(df):
            
    sns.boxplot(x='LCD Class', y='Percent Identity', data=df, showfliers=False, color='#1f77b4')
    sns.stripplot(x='LCD Class', y='Percent Identity', data=df, color='0.2', jitter=True, alpha=0.5, s=2)
    plt.plot([-1, 20], [80, 80], linestyle='--', color='0.5')
    plt.xticks(fontname='Arial', fontsize=16)
    plt.yticks(fontname='Arial', fontsize=16)
    plt.xlabel('LCD Class', fontname='Arial', fontsize=18)
    plt.ylabel('Percent Sequence Identity', fontname='Arial', fontsize=18)
    plt.xlim(-0.55, 18.5)
    fig = plt.gcf()
    fig.set_size_inches(8, 5)
    plt.savefig('Fig S15A - Main_Yeast_LCDs_PercentIdentity_Boxplot.tiff', bbox_inches='tight', dpi=600)
    plt.close()
    

if __name__ == '__main__':
    main()