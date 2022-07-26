
"""
Description:
    LCD-Composer is a composition-based method for identifying low-complexity domains in protein sequences.
    Refer to Cascarina et. al. (2021) (https://academic.oup.com/nargab/article/3/2/lqab048/6285187) for a complete 
    description of the original LCD-Composer algorithm, which is used in original and modified form here.
    
    This script is a modification/augmentation of the original LCD-Composer algorithm, developed as 
    part of the LCD-Composer webserver.
    
    For a description of modifications and added functionality, refer to Cascarina and Ross (2022) 
    (LINK TO BIOINFORMATICS WHEN PUBLISHED) and the LCD-Composer webserver HELP page at 
    http://lcd-composer.bmb.colostate.edu:9000/help.
    Documentation for LCD-Composer is available at https://github.com/RossLabCSU/LCD-Composer.
============================================================================================
License_info:
    LCD-Composer is subject to the terms of the GPLv3.0 license. For a complete description of 
    license terms, please see the license at https://github.com/RossLabCSU/LCD-Composer.
"""

__author__ = 'Sean M Cascarina'
__copyright__ = 'Copyright 2022'
__credits__ = ['Sean M Cascarina']
__license__ = 'GPLv3.0'
__version__ = '1.0'
__maintainer__ = 'Sean M Cascarina'
__email__ = 'Sean.Cascarina@colostate.edu'


all_aas = 'ACDEFGHIKLMNPQRSTVWY'


def analyze_cooccurrence(jobfile_ids, output):
    
    df, lcd_type_df = retrieve_jobfile_data(jobfile_ids)
    
    for prot in df:
        if len(set(df[prot]['Jobfiles'])) < len(jobfile_ids):   # SKIPS PROTEINS THAT DID NOT HAVE AT LEAST ONE LCD IN EACH jobfile_id
            continue
            
        overlap_checks = crosscheck_LCD_bounds(df[prot]['LCD Bounds'])
        lcd_types = [lcd_type_df[jobfile] for jobfile in df[prot]['Jobfiles']]
        
        for i, line_data in enumerate(df[prot]['Lines']):
            line_data += [lcd_types[i], overlap_checks[i]]
            output.write('\t'.join(line_data) + '\n')
            
    return output
    

def retrieve_jobfile_data(jobfile_ids):

    df = {}
    lcd_type_df = {}
    for jobfile_id in jobfile_ids:
        filename = jobfile_id + '_LCDcomposer_RESULTS.tsv'
        h = open(filename)
        
        # READ METADATA
        line = '>DummyString'
        while line.startswith('>'):
            line = h.readline()
            if line.startswith('>Amino Acid'):
                junk, aas = line.rstrip().split(': ')
                aas = aas.replace('[', '').replace(']', '')
                aas = aas[:-1].replace("'", '').split(', ')
                aa_string = '_'.join(aas)
                lcd_type_df[jobfile_id] = aa_string
        
        header = h.readline().rstrip().split('\t')
        
        for line in h:
            items = line.rstrip().split('\t')
            
            # HANDLES DATA FROM OPTION2 OUTPUT
            if header[0] == 'LCD Similarity Rank':
                items = items[2:]
                
            uniprot = items[1]
            start, end = items[3][1:-1].split('-')
            start, end = int(start), int(end)

            df[uniprot] = df.get(uniprot, {'Lines':[], 'LCD Bounds':[], 'Jobfiles':[]})
            df[uniprot]['Lines'].append(items[:6])
            df[uniprot]['LCD Bounds'].append( (start, end) )
            df[uniprot]['Jobfiles'].append( jobfile_id )
            
        h.close()
        
    return df, lcd_type_df

        
def crosscheck_LCD_bounds(bounds):

    overlap_checks = []
    for ind1, bound1 in enumerate(bounds):
        overlap = 'No'
        for ind2, bound2 in enumerate(bounds):
            if ind1 == ind2:
                continue
            
            if (bound1[0] >= bound2[0] and bound1[0] <= bound2[1]) or (bound1[1] >= bound2[0] and bound1[1] <= bound2[1]) or (bound1[0] >= bound2[0] and bound1[1] <= bound2[1]) or (bound1[0] <= bound2[0] and bound1[1] >= bound2[1]):
                overlap = 'Yes'
                
        overlap_checks.append(overlap)
        
    return overlap_checks
    

def get_params(args):
    """Gather, define, and validate user parameters."""
    
    jobfile_ids = args.jobfile_ids
    allowed_fileID_chars = '1234567890_'
    if sum([1 for char in jobfile_ids if char in allowed_fileID_chars]) != len(jobfile_ids):
        print('\nError:\nYour Job ID numbers can only contain numbers 0-9 and underscores.\n')
        exit()
        jobfile_ids = {''}
    else:
        jobfile_ids = set( jobfile_ids.split('_') )
        if len(jobfile_ids) < 2:
            print('\nError:\nThere appears to be only one Job ID number submitted. Please check the formatting of your submitted Job IDs (one per line, no spaces) and ensure that you are submitting 2-3 Job IDs.\n')

    failed_fileopen_test = False
    is_autodetect = False
    is_option3 = False
    for jobfile_id in jobfile_ids:
        filename = jobfile_id + '_LCDcomposer_RESULTS.tsv'
        try:
            h = open(filename)
            for line in h:
                if 'Auto-detect' in line:
                    is_autodetect = True
                if line.startswith('Protein ID'):
                    is_option3 = True
            h.close()
        except:
            failed_fileopen_test = True
            
    if failed_fileopen_test:
        print('\nError:\nOne or more of your LCD search results files could not be found. Please ensure that the file IDs are correct, that the file name has not been altered, and that the files exist in the same folder/directory as this script.\n')
        exit()
        
    if is_autodetect:
        print('\nError:\nOption 5 is not suitable for comparing results files that are generated using "Auto-detect" mode. If you are interested in exploring what other types of LCDs are in your set of LCD-containing proteins (without a specific prior hypothesis), we recommend using "Auto-detect mode" in Option 4 with your proteins of interest. This will search for each of the 20 primary classes of LCDs (one for each amino acid) within your proteins of interest.\n')
        exit()
        
    if is_option3:
        print('\nError:\nOption 5 is not suitable for comparing results files that are generated using Option 3 (single-protein analysis) since Option 3 itself can easily be used to discover multiple types of LCDs in your protein of interest.\n')
        exit()
    
    return jobfile_ids
    
    
def output_params_header(output, random_id, jobfile_ids):
    
    output.write('>*RUNTIME PARAMETERS*\n')
    output.write('>Job ID of Co-Occurrence Results File: ' + str(random_id) + '\n')
    output.write('>IDs of Files Compared: ' + ', '.join(jobfile_ids) + '\n')
    
    return output

    
def main(args):

    jobfile_ids = get_params(args)
    random_id = random.randint(0, 10000000)

    results_file = str(random_id) + '_LCD_Co-Occurrence_RESULTS.tsv'
    output = open(results_file, 'w')
    output_params_header(output, random_id, jobfile_ids)
    output.write('\n')
    output.write('\t'.join( ['Protein Description','UniProt ID','Domain Sequence','Domain Boundaries','Final Domain Composition','Final Domain Linear Dispersion','Amino Acid(s) Used in Original LCD Search','LCD Overlap with LCD from Other Results File?']) + '\n')
    
    output = analyze_cooccurrence(jobfile_ids, output)
    
    output.close()
        

def get_args(arguments):
    parser = argparse.ArgumentParser(description='Comparison of LCD search results files to identify multi-LCD proteins.')
    
    parser.add_argument('jobfile_ids', help="""Your LCD-Composer results file IDs, underscore ("_") separated.""")

    args = parser.parse_args(arguments)
    
    return args
    

if __name__ == '__main__':
    import sys, argparse, datetime, math, statistics, random, os
    args = get_args(sys.argv[1:])
    main(args)