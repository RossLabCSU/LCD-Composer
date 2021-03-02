
amino_acids = 'ACDEFGHIKLMNPQRSTVY' #W OMITTED BECAUSE NO W-RICH LCDs ARE IDENTIFIED BY THE STANDARD 40% COMPOSITION THRESHOLD

def main():

    h = open('SecondaryAA_LCD_Statistics_Uncorrected_Pvals.csv')
    output = open('TableS10_SecondaryAA_LCD_Statistics_Corrected_Pvals.csv', 'w')
    explanation_line = h.readline()
    header = h.readline().rstrip().split(',')
    header = header[:3] + ['Bonferroni-corrected P-value'] + header[3:]
    output.write(','.join(header) + '\n')
    
    for aa in amino_acids:
        pvals = []
        master_items = []
        
        #CALCULATE NUMBER OF NANS AND STORE EACH LINE
        num_nans = 0
        for i in range(19):
            items = h.readline().rstrip().split(',')
            master_items.append(items)
            pval = items[2]
            if pval == 'nan':
                num_nans += 1
                
        #LOOP THROUGH LINES, CALCULATE BONF CORRECTION, AND OUTPUT NEW LINES
        for items in master_items:
            pval = items[2]
            if pval != 'nan':
                pval = float(pval)
                output.write(','.join( items[:3] + [str( min(pval*(19-num_nans), 1) )] + items[3:] ) + '\n')
            else:
                output.write(','.join( items[:3] + ['nan'] + items[3:] ) + '\n')
                
    h.close()
    output.close()
    
if __name__ == '__main__':
    main()