
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
import numpy as np
import statistics
import matplotlib.pyplot as plt
import seaborn as sns
import scipy as sp
import pandas as pd


def get_data_v2():
    import datetime

    df = {'Computation Time':[],
        'Window Size':[],
        'Composition Threshold':[],
        'Linear Dispersion Threshold':[],
        'Category':[],
        'Amino Acid':[]}
    
    cats = ['Composition Threshold', 'Linear Dispersion Threshold', 'Window Size']
    index = 0
    all_times = []

    for file in ['ComputationTimes_CompThresholdVaried.txt', 'ComputationTimes_DispThresholdVaried.txt', 'ComputationTimes_WindowSizeVaried.txt']:

        h = open(file)

        win_sizes = []
        comp_thresholds = []
        disp_threshold = []
        aas = []
        win_size = 20
        for line in h:
            if 'LCD-Composer' in line:
                scer, aa, comp, winsize, disp, lcdcomposer, results = line.rstrip().split('_')
                comp_threshold, temp = comp.split('-')
                comp_threshold = int(comp_threshold)
                
                #SKIP EXCESSIVELY LOW COMPOSITION THRESHOLDS
                if comp_threshold in [10, 20]:
                    h.readline()
                    h.readline()
                    h.readline()
                    continue
                    
                win_size, temp = winsize.split('-')
                win_size = int(win_size)
                disp_threshold, temp = disp.split('-')

            if 'Start time' in line:
                items = line.rstrip().split(' ')
                start_dt = items[2] + ' ' + items[3]
                start_dtobj = datetime.datetime.strptime(start_dt, '%Y-%m-%d %H:%M:%S.%f')
                
                #PARSE NEXT LINE, WHICH IS END TIME
                items = h.readline().rstrip().split(' ')
                end_dt = items[2] + ' ' + items[3]
                end_dtobj = datetime.datetime.strptime(end_dt, '%Y-%m-%d %H:%M:%S.%f')
                
                delta_t = end_dtobj - start_dtobj
                seconds = delta_t.total_seconds()

                df['Computation Time'].append(seconds)
                df['Window Size'].append(win_size)
                df['Composition Threshold'].append(comp_threshold)
                df['Linear Dispersion Threshold'].append(disp_threshold)
                df['Category'].append(cats[index])
                df['Amino Acid'].append( aa )

                all_times.append(seconds)
        index += 1
        h.close()
            
    new_df = {'Computation Time':[],
        'Category':[],
        'Standard Deviation':[],
        'Amino Acid':[]}

    data = pd.DataFrame.from_dict(df)
    
    for cat in cats:
        cat_df = data[data['Category'] == cat]
        for aa in amino_acids:
            temp_df = cat_df[cat_df['Amino Acid'] == aa]
            times = list(temp_df['Computation Time'])
            if len(times) > 0:
                new_df['Computation Time'].append( sp.nanmean(times) )
            else:
                new_df['Computation Time'].append( np.nan )
            new_df['Standard Deviation'].append( np.nanstd(times) )
            new_df['Category'].append(cat)
            new_df['Amino Acid'].append(aa)
            
    new_df = pd.DataFrame.from_dict(new_df)
        
    grouped_barplot(new_df, cats)


def grouped_barplot(df, cats):

    palette = sns.color_palette('colorblind')
    xvals = [x for x in range(len(amino_acids))]
    offsets = [-0.25, 0, 0.25]
    
    for i in range(len(cats)):
        cat = cats[i]
        small_df = df[df['Category'] == cat]
        plt.bar([x+offsets[i] for x in xvals], small_df['Computation Time'].values, width=0.25, label="{} {}".format('Category', cat), yerr=small_df['Standard Deviation'].values, color=palette[i])
        
    plt.xlabel('Amino Acid', fontname='Arial', fontsize=18)
    plt.ylabel('Computation Time\n(seconds)', fontname='Arial', fontsize=18)
    plt.xticks(xvals, list(amino_acids), fontname='Arial', fontsize=16)
    plt.yticks(fontname='Arial', fontsize=16)
    
    max_val = np.nanmax( [df['Computation Time'][i] + df['Standard Deviation'][i] for i in range(len(df['Computation Time']))] )
    plt.ylim(0, max_val+2)
    plt.xlim(-0.5, 19.5)
    
    leg = plt.legend(labels=['Composition Threshold', 'Linear Dispersion Threshold', 'Window Size'], prop={'family':'Arial', 'size':12}, loc=2, handletextpad=0.3, borderpad=0.4, handlelength=2, labelspacing=0.15)
    leg.set_title('Varied Parameter', prop={'family':'Arial', 'size':14})
    
    fig = plt.gcf()
    fig.set_size_inches(10, 5)
    plt.savefig('Fig S6A - Computation Times_Yeast_LCD-Composer_ParametersVaried.tiff', bbox_inches='tight', dpi=600)
    plt.close()

    
if __name__ == '__main__':
    get_data_v2()