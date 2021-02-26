
import matplotlib.pyplot as plt
import seaborn as sns
import statistics
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
import pandas as pd
import numpy as np
aa_names = {'A':'Alanine (A)', 'C':'Cysteine (C)', 'D':'Aspartic Acid (D)', 'E':'Glutamic Acid (E)',
        'F':'Phenylalanine (F)', 'G':'Glycine (G)', 'H':'Histidine (H)', 'I':'Isoleucine (I)', 
        'K':'Lysine (K)', 'L':'Leucine (L)', 'M':'Methionine (M)', 'N':'Asparagine (N)',
        'P':'Proline (P)', 'Q':'Glutamine (Q)', 'R':'Arginine (R)', 'S':'Serine (S)',
        'T':'Threonine (T)', 'V':'Valine (V)', 'W':'Tryptophan (W)','Y':'Tyrosine (Y)'}

def main():
    
    total_domains_wideform = {}
    domain_lens_wideform = {}

    domain_lengths = {'Composition Threshold':[],
                    'Window Size':[],
                    'Dispersion Threshold':[],
                    'Length':[],
                    'Amino Acid':[],
                    'Varied Parameter':[]}
                            
    # total_domains_longform, domain_lengths, total_domains_wideform, domain_lens_wideform = get_my_data(total_domains_longform, domain_lengths, total_domains_wideform, domain_lens_wideform)
    # total_domains_longform, domain_lengths, total_domains_wideform = get_my_data(domain_lengths, total_domains_wideform)
    # domain_lengths = get_my_data(domain_lengths)

    # domain_lengths = pd.DataFrame.from_dict(domain_lengths)

    # print(domain_lengths)

    for category in ['Window Size', 'Dispersion Threshold', 'Composition Threshold']:
        domain_lengths = {'Composition Threshold':[],
                        'Window Size':[],
                        'Dispersion Threshold':[],
                        'Length':[],
                        'Amino Acid':[],
                        'Varied Parameter':[]}
        domain_lengths = get_my_data(domain_lengths, category)

        # domain_lengths = pd.DataFrame.from_dict(domain_lengths)
        print(domain_lengths)
        if category == 'Window Size':
            len_df = pd.DataFrame.from_dict(domain_lengths)
            len_df = len_df.loc[(len_df['Dispersion Threshold']=='0.5') & (len_df['Composition Threshold']==40)]
        elif category == 'Dispersion Threshold':
            len_df = pd.DataFrame.from_dict(domain_lengths)
            len_df = len_df.loc[(len_df['Window Size']==20) & (len_df['Composition Threshold']==40)]
        else:
            len_df = pd.DataFrame.from_dict(domain_lengths)
            len_df = len_df.loc[(len_df['Window Size']==20) & (len_df['Dispersion Threshold']=='0.5')]
        
        huesplit_barcharts(len_df, category)

    
    
# def get_my_data(total_domains, domain_lengths, totals_wideform):
def get_my_data(domain_lengths, category):

    # for aa in all_aas:
        # amino_acids = [aa]
        # combined_aas = aa
        # other_aas = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ*_'.replace(aa, '')
        # for win_size in range(10, 110, 10):
            # for comp_threshold in range(10, 110, 10):
                # comp_thresholds = [comp_threshold]
                # ign_disp = comp_threshold + ((100 - comp_threshold)/2)
                # for disp_threshold in [x/10 for x in range(0, 11)]:
                
    # if category == 'Window Size':
        # files = ['Scerevisiae_' + aa + '_' + str(comp_threshold) + '-CompositionThreshold_' + str(win_size) + '-WindowSize_' + str(disp_threshold).replace('.', 'p') + '-LinearDispersionThrshold_LCD-Composer_RESULTS.tsv')
        # comp_threshold = '40'
        # disp_threshold = '0p5'
        # varied_params = [
    # elif 
                
    for aa in amino_acids:
        if category == 'Window Size':
            files = ['Scerevisiae_'+aa+'_40-CompositionThreshold_'+str(win_size)+'-WindowSize_0p5-LinearDispersionThreshold_LCD-Composer_RESULTS.tsv' for win_size in range(10, 110, 10)]
            comp_thresholds = [40 for i in range(len(files))]
            disp_thresholds = ['0p5' for i in range(len(files))]
            win_sizes = [x for x in range(10, 110, 10)]
        elif category == 'Dispersion Threshold':
            files = ['Scerevisiae_'+aa+'_40-CompositionThreshold_20-WindowSize_'+disp_threshold+'-LinearDispersionThreshold_LCD-Composer_RESULTS.tsv' for disp_threshold in ['0p0', '0p1', '0p2', '0p3', '0p4', '0p5', '0p6', '0p7', '0p8', '0p9', '1p0']]
            comp_thresholds = [40 for i in range(len(files))]
            disp_thresholds = [x for x in ['0p0', '0p1', '0p2', '0p3', '0p4', '0p5', '0p6', '0p7', '0p8', '0p9', '1p0']]
            win_sizes = [20 for i in range(len(files))]
        else:
            files = ['Scerevisiae_'+aa+'_'+str(comp_threshold)+'-CompositionThreshold_20-WindowSize_0p5-LinearDispersionThreshold_LCD-Composer_RESULTS.tsv' for comp_threshold in range(10, 110, 10)]
            comp_thresholds = [x for x in range(10, 110, 10)]
            disp_thresholds = ['0p5' for i in range(len(files))]
            win_sizes = [20 for i in range(len(files))]
                

        for i in range(len(files)):
        # for file in files:
    # for aa in amino_acids:
        # for win_size in range(10, 110, 10:
        # totals_wideform['Amino Acid'].append(aa)
            # for comp_threshold in range(10, 110, 10):
                # for disp_threshold in [x/10 for x in range(0, 11)]:

                    # h = open('uniprot-proteome_UP000002311_YeastProteome_' + aa + '_20aaWindow_' + str(comp_threshold) + 'CompThreshold_0,5DispThreshold_RESULTS.csv')
            # h = open('Scerevisiae_' + aa + '_' + str(comp_threshold) + '-CompositionThreshold_' + str(win_size) + '-WindowSize_' + str(disp_threshold).replace('.', 'p') + '-LinearDispersionThrshold_LCD-Composer_RESULTS.tsv')
            
            # NEED TO REMOVE THIS TRY/EXCEPT AFTER TESTING AND ONLY HAVE THE h = open(files[i]) LINE==================================================================================================================================================
            # try:
                # h = open(files[i])
            # except:
                # h = open(files[i] + '.tsv')
            h = open(files[i])
            #======================================================================================================================================================================================================================================
            
                
            # h = open('uniprot-proteome_UP000002311_YeastProteome_' + aa + '_' + str(win_size) + 'aaWindow_40CompThreshold_0,3SepThreshold_Results_WithMustContain.csv')
            for j in range(7):
                h.readline()
            header = h.readline()
            counts = 0
            lengths = []
            for line in h:
                id, seq, boundaries, final_comp, final_stddev, *remainder = line.rstrip().split('\t')
                counts += 1
                lengths.append(len(seq))
                
            h.close()

            domain_lengths['Composition Threshold'].append( comp_thresholds[i]  )
            domain_lengths['Window Size'].append( win_sizes[i]  )
            domain_lengths['Dispersion Threshold'].append( disp_thresholds[i].replace('p', '.')  )
            domain_lengths['Amino Acid'].append( aa )
            domain_lengths['Varied Parameter'].append( category )
            if len(lengths) > 0:
                domain_lengths['Length'].append( statistics.mean(lengths) )
            else:
                domain_lengths['Length'].append( 0 )
                
            h.close()
    # print(totals_wideform)
                
    return domain_lengths
    
    
def huesplit_barcharts(domain_lengths, category):
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf', '0.8']
    
    # sns.barplot(x=domain_lengths['Amino Acid'], y=domain_lengths['Length'], hue=domain_lengths['SepThreshold'], palette=colors)
    # sns.barplot(x=domain_lengths['Amino Acid'], y=domain_lengths['Length'], hue=domain_lengths['Window Size'], palette=colors)
    print(domain_lengths)
    # sns.barplot(x=domain_lengths['Amino Acid'], y=domain_lengths['Length'], hue=domain_lengths['Composition Threshold'], palette=colors)
    sns.barplot(x=domain_lengths['Amino Acid'], y=domain_lengths['Length'], hue=domain_lengths[category], palette=colors)
    
    plt.xticks(fontname='Arial', fontsize=18)
    plt.yticks(fontname='Arial', fontsize=18)
    plt.xlabel('Amino Acid', fontname='Arial', fontsize=20)
    plt.ylabel('Average Domain Length', fontname='Arial', fontsize=20)
    leg = plt.legend(prop={'family':'Arial', 'size':14}, loc=2, bbox_to_anchor=(1.0, 1.02))
    # leg.set_title('Linear\nDispersion\nThreshold', prop={'family':'Arial', 'size':14})
    leg.set_title(category.replace(' ', '\n'), prop={'family':'Arial', 'size':14})
    plt.tight_layout()
    fig = plt.gcf()
    fig.set_size_inches(12, 6)
    # plt.savefig('AveDomainLengths_FixedDispersion0,5_CompThresholdVaried.tiff', bbox_inches='tight', dpi=600)
    if category == 'Window Size':
        plt.savefig('Fig S5A - AveDomainLengths_' + category + ' Varied.tiff', bbox_inches='tight', dpi=600)
    elif category == 'Dispersion Threshold':
        plt.savefig('Fig S5B - AveDomainLengths_' + category + ' Varied.tiff', bbox_inches='tight', dpi=600)
    else:
        plt.savefig('Fig S5C - AveDomainLengths_' + category + ' Varied.tiff', bbox_inches='tight', dpi=600)
    plt.close()
    # plt.show()
    
def overlaid_barcharts(domain_lengths):

    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf', '0.5']
    
    index = 0
    for sep_threshold in [str(round(x/10, 1)) for x in range(0, 11)]:
    
        xvals = []
        yvals = []
        labels = []
        for i in range(len(amino_acids)):
            xvals.append(i)
            labels.append(amino_acids[i])
            small_df = domain_lengths[domain_lengths['Amino Acid'] == amino_acids[i]]
            # print(small_df)
            # small_df = small_df[small_df['SepThreshold'] == sep_threshold.replace('.', ',')]
            small_df = small_df[small_df['SepThreshold'] == sep_threshold]
            # print(small_df)
            yval = list(small_df['Length'])[0]
            # yvals.append( statistics.mean( small_df[small_df['SepThreshold'] == sep_threshold] ) )
            yvals.append( yval )
            
    
        # small_df = domain_lengths[domain_lengths['SepThreshold'] == sep_threshold]
        # sns.barplot(x=small_df['Amino Acid'], y=[statistics.mean(small_df[
        sns.barplot(x=xvals, y=yvals, color=colors[index])
        
        index += 1
    plt.xticks([x for x in range(20)], labels=labels)
        
    plt.show()
        
    
def plotting(df):

    print(df)
        # plt.scatter(new_df['Composition'], new_df['Spacing'])
        # plt.show()
        
    aa_index = 1
    #fig = plt.figure(figsize=(8,11))
    #fig, ((ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8), (ax9, ax10, ax11, ax12), (ax13, ax14, ax15, ax16), (ax17, ax18, ax19, ax20)) = plt.subplots(5,4,sharex='col', sharey='row', figsize=(8,11))
    fig, ax = plt.subplots(5,4,sharex=True, sharey=False, figsize=(8,10))
    # lgnd = plt.legend(['Positives', 'Ambiguous', 'Negatives'], loc=1)
    big_subplot = fig.add_subplot(111)
    # big_subplot.set_xlabel('Category', fontname='Arial')
    # big_subplot.set_ylabel('Percent Composition', fontname='Arial')
    big_subplot.spines['top'].set_color('none')
    big_subplot.spines['bottom'].set_color('none')
    big_subplot.spines['left'].set_color('none')
    big_subplot.spines['right'].set_color('none')
    big_subplot.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)
    # big_subplot.legend(['Positives', 'Negatives'], loc=1)
    
    for aa in amino_acids:

        new_df = df[df['Amino Acid'] == aa]
        # new_df = new_df.to_dict()
        
        # small_df = { 'Percent Composition' : [],
                    # 'Dataset' : [],
                    # 'Category' : []}
        # pos_comps = []
        # neg_comps = []
        # index = 0
        # for seq in df['Sequence']:
            # perc_comp = seq.count(aa) / len(seq) * 100
            # small_df['Category'].append(df['Category'][index])
            # small_df['Percent Composition'].append(perc_comp)
            # small_df['Dataset'].append(df['Dataset'][index])
            # index += 1
            
        #PLOTTING
        ax = plt.subplot(5,4, aa_index)
        # plt.scatter(new_df['Composition'], new_df['Spacing'], s=4, alpha=0.4)
        print(list(new_df['SepThreshold']))
        print(list(new_df['Counts']))
        # plt.bar(new_df['SepThreshold'], height=new_df['Counts'])
        sns.barplot(x=new_df['SepThreshold'], y=new_df['Counts'], color='#1f77b4')
        # plt.bar(new_df['Counts'])
        # ax.legend.remove()
        plt.ylabel('')
        plt.xlabel('')

        plt.title(aa_names[aa], fontsize=12, fontname='Arial')
        plt.xlim(-1, 10)
        # plt.ylim(-0.05, 1.05)
        # plt.xticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0], fontname='Arial', fontsize=14)
        # plt.xticks([2, 4, 6, 8, 10], fontname='Arial', fontsize=14, rotation=90)
        plt.xticks([1, 3, 5, 7, 9], fontname='Arial', fontsize=14, rotation=90)
        
        if aa == 'W':
            plt.yticks([0, 1], labels=['0', '1'], fontname='Arial', fontsize=14)
        else:
            plt.yticks(fontname='Arial', fontsize=14)
        
        # plt.tight_layout()
        #ax.set_xticklabels([])
        # plt.tick_params(
            # axis='x',          # changes apply to the x-axis
            # which='both',      # both major and minor ticks are affected
            # bottom=False,      # ticks along the bottom edge are off
            # top=False,         # ticks along the top edge are off
            # labelbottom=False) # labels along the bottom edge are off
        #plt.grid(color='0.95')
        #plt.axis([-0.5, 2.5])
        #plt.xticks(np.arange(0,105,10))
        #plt.yticks(np.arange(0,105,10))
        
        #set labels for axes on outside subplots only
        # if aa_index in (2,3,4,6,7,8,10,11,12,14,15,16,18,19,20):
            # ax.set_yticklabels([])
        # if aa_index in (1,5,9,13,17):
            # ax.set_yticklabels(['0.0','0.2','0.4','0.6','0.8','1.0'])
        if aa_index in range(1,17):
            ax.set_xticklabels([])
        if aa_index in range(17,21):
            ax.set_xticklabels(['0.2','0.4','0.6','0.8','1.0'])
        #plt.show()
        aa_index += 1

    #fig.text(0.5, 0.0, 'Category', ha='center', fontname='Arial')
    fig.text(-0.03, 0.5, 'Number of LCDs', va='center', rotation='vertical', fontname='Arial', fontsize=16)
    fig.text(0.5, -0.02, 'Linear Dispersion Threshold', ha='center', fontname='Arial', fontsize=16)
    plt.tight_layout(pad=0.2)
    # plt.savefig('TotalNumLCDs_20aaWindow_0,3Composition_SUBPLOTS.tiff', bbox_inches ='tight', dpi=600)
    plt.close()  
    
    
def heatmap_total_nums(df):
    
    df.set_index('Amino Acid', inplace=True)
    for key in df.columns:
        print(key, repr(key), type(key))
    print(df.columns)
    # df.sort_values(by=['1.0'], axis=1, inplace=True, ascending=False)
    df.sort_values('0.0', inplace=True, ascending=False)
    
    log_df = np.log10(df.replace(0, np.nan))
    log_df.replace(np.nan, 0, inplace=True)
    
    colors = sns.color_palette("coolwarm", 1000)
    # ax = sns.heatmap(df, cmap=colors, annot=True, fmt='#', annot_kws={'family':'Arial', 'size':12}, cbar_kws={'label':'# of LCDs'})
    
    # ax = sns.heatmap(df, cmap=colors, annot=True, fmt='#', annot_kws={'family':'Arial', 'size':12})
    ax = sns.heatmap(log_df, cmap=colors, annot=df, fmt='#', annot_kws={'family':'Arial', 'size':12})
    
    # plt.xticks([x for x in range(1, 11)], labels = [str(round(x/10, 1)) for x in range(1,11)], fontname='Arial', fontsize=16)
    plt.xticks(fontname='Arial', fontsize=16)
    plt.yticks(fontname='Arial', fontsize=16, rotation=0)
    plt.xlabel('Linear Dispersion Threshold', fontname='Arial', fontsize=18)
    plt.ylabel('Amino Acid', fontname='Arial', fontsize=18)
    
    #COLORBAR
    # https://stackoverflow.com/questions/37233108/seaborn-change-font-size-of-the-colorbar
    # https://stackoverflow.com/questions/34820239/seaborn-heatmap-colorbar-label-as-percentage
    cbar = ax.collections[0].colorbar
    labels = cbar.ax.yaxis.get_ticklabels()
    # cbar.set_ticks([0, .2, .75, 1])
    # cbar.ax.set_title('# of LCDs', fontname='Arial', fontsize=16, ha='left')
    cbar.ax.set_title('log( # of LCDs )', fontname='Arial', fontsize=16, ha='center')
    cbar.ax.set_yticklabels(labels, fontname='Arial', fontsize=16)
    # cbar.ax.yaxis.label.set_size(20)
    
    # ax = plt.gca()
    # im = ax.imshow(df)
    # clb = plt.colorbar(ax=ax)
    # clb.ax.set_title('# of LCDs')
    # labels = clb.ax.yaxis.get_ticklabels()
    # clb.ax.set_yticks(labels, fontname='Arial', fontsize=16)
    
    # plt.show()
    fig = plt.gcf()
    fig.set_size_inches(8, 8)
    # plt.savefig('Total Number of LCDs_SepThresholdScanning_InverseDispersion_WithMustContain_HEATMAP.tiff', bbox_inches='tight', dpi=600)
    plt.close()
    
    
def plot_total_nums(total_domains, thresholds):
    
    totals = []
    xvals = [x for x in range(1, len(totals)+1)]
    for key in sorted(list(total_domains.keys())):
        totals.append(total_domains[key])
        print(key, total_domains[key])
        
    sns.barplot(thresholds, totals)
    # plt.xticks(labels=[str(round(x, 1)) for x in thresholds])
    plt.show()
    
    
def plot_length_dist(domain_lengths):
    
    sns.distplot(domain_lengths, hue='SepThreshold', bins=30)
    plt.show()
    
def plot_length_dist_shortform(domain_lens):
    
    for key in sorted(list(domain_lens.keys())):
        sns.distplot(domain_lens[key], bins=30)
    
    plt.show()

if __name__ == '__main__':
    main()