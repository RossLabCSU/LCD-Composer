
gaf_df = {'UP000006548_3702':'goa_arabidopsis.gaf', 'UP000009136_9913':'goa_cow.gaf', 'UP000001940_6239':'goa_worm.gaf', 'UP000002254_9615':'goa_dog.gaf', 'UP000000437_7955':'goa_zebrafish.gaf', 'UP000002195_44689':'goa_dicty.gaf', 'UP000000803_7227':'goa_fly.gaf', 'UP000000539_9031':'goa_chicken.gaf', 'UP000005640_9606':'goa_human.gaf', 'UP000000589_10090':'goa_mouse.gaf', 'UP000002311_559292':'goa_yeast.gaf', 'UP000008227_9823':'goa_pig.gaf', 'UP000002494_10116':'goa_rat.gaf'}

def main(args):
    
    jobfile_id, lcd_resultsfile, fasta_file, gaf_file = get_params(args)
    
    lcd_prots = get_LCDprots(jobfile_id)

    go_obo = 'go-basic.obo'
    go = obo_parser.GODag(go_obo, prt=sys.stdout)
    
    # RUN CHECK TO ENSURE THAT THE GAF FILE MATCHES THE PROTEOME FILE IF THE USER IS EXAMINING A SUPPORTED ORGANISM==========
    extensions = ['.fsa', '.fasta', '.fnn', '.FASTA', '.fa']
    fasta_str = fasta_file[:]
    if fasta_str.startswith('.\\'):
        fasta_str = fasta_str[2:]
    check_ext = [1 if ext in fasta_str else 0 for ext in extensions]
    extension = extensions[ check_ext.index(1) ]
    fasta_str = fasta_str.replace(extension, '')
    
    if fasta_str in gaf_df and not (gaf_file != gaf_df[fasta_str] or gaf_file != '.\\' + gaf_df[fasta_str]):
        print('\n#==========================================#')
        print('WARNING: Based on our internal database, your specified GAF file does not appear to match the proteome file used to define the background. Below is a list of expected filenames:')
        for key, value in gaf_df.items():
            print('Proteome file: ' + key + '\t\tCorresponding GAF file: ' + value)
        print('\nIt is possible that you are using an alternative GAF file and/or proteome. Please double-check that these files correspond to the same organism.')
        print('#==========================================#')
    # =================================

    background_prots = get_background_prots( fasta_file )
    assoc = get_assoc(gaf_file)

    methods = ["bonferroni", "sidak", "holm", "fdr"]
    g = GOEnrichmentStudy(background_prots, assoc, go, propagate_counts=False, alpha=0.05)
    g_res = g.run_study(lcd_prots)
        
    output_GOresults(jobfile_id, fasta_str, g_res)

    
def get_LCDprots(jobfile_id):

    lcd_prots = set()
    filename = jobfile_id + '_LCDcomposer_RESULTS.tsv'
    try:
        h = open(filename)
    except:
        print('\nError:\nYour LCD search results file could not be found. Please ensure that the file ID is correct, that the file name has not been altered, and that the file exists in the same folder/directory as this script.\n')
        exit()
    
    # READ METADATA
    line = '>DummyString'
    while line.startswith('>'):
        line = h.readline()

    header = h.readline().rstrip().split('\t')
    
    for line in h:
        items = line.rstrip().split('\t')
        
        # HANDLES DATA FROM OPTION2 OUTPUT
        if header[0] == 'LCD Similarity Rank':
            items = items[2:]
            
        uniprot = items[1]
        lcd_prots.add(uniprot)

    h.close()
        
    return lcd_prots

        
def get_background_prots(proteome):
    
    prots = []
    h = open(proteome)
    for id, seq in fasta_parser(h):
        try:
            junk, uniprot_id, junk = id.split('|')
        except:
            print('Error:\nYour proteome file does not appear to be a UniProt proteome. This program only supports select UniProt proteomes.\n')
            exit()
        prots.append( uniprot_id )
        
    return prots
    
    
def get_assoc(abbrev):
    """
    Reads the proper gene association file (GAF).
    Return is a dictionary with Protein IDs as keys, and a set of associated GO terms as values (formal Python set)
    """

    h = open(abbrev)

    assoc = {}

    #POPULATE 'assoc' DICTIONARY TO PASS TO GOEnrichmentStudy() MODULE IN GOATOOLS
    for line in h:
        if line.startswith('!'):
            continue
        items = line.rstrip().split('\t')
        uniprot_id = items[1]
        assoc[uniprot_id] = assoc.get(uniprot_id, set())
        assoc[uniprot_id].add(items[4])

    return assoc
    
    
def fasta_parser(file):
    """Parses each instance in a FASTA formatted file into a gene id and a sequence.
    
    Yields:
        id, seq (both strings)
    """
    #INITIALIZES GENE AND SEQ FOR FIRST ITERATION
    gene, seq = '', []
    
    #LOOPING THROUGH FILE USING GENERATOR
    for line in file:
        line = line.rstrip()
        if len(line) == 0:
            continue
            
        #YIELDS GENE AND SEQ IF THE NEXT ID IS REACHED
        if line.startswith('>'):
            if gene != '':
                yield (gene[1:], ''.join(seq))
            gene, seq = line, []
        else:
            seq.append(line)
            
    #YIELDS FINAL INSTANCE IN FASTA FILE
    if gene != '':
        yield (gene[1:], ''.join(seq))
    
    
def get_params(args):

    jobfile_id = args.jobfile_id
    lcd_resultsfile = jobfile_id + '_LCDcomposer_RESULTS.tsv'
    fasta_file = args.fasta_file
    gaf_file = args.gaf_file
    
    for file in [lcd_resultsfile, fasta_file, gaf_file]:
        try:
            h = open(file)
            h.close()
        except:
            print('\nError:\nYour file "' + file + '" could not be found. Please ensure that it is correctly named and is in the same folder/directory as this script. You must also include the correct file extension in the file names.\n') 
            exit()
            
    return jobfile_id, lcd_resultsfile, fasta_file, gaf_file
    
    
def crosscheck_proteome_file(lcd_prots, background_prots):

    valid_prots = [prot for prot in lcd_prots if prot in background_prots]
    if len(valid_prots) == 0:
        print('\nError:\nNone of the protein IDs from your job file could be found in the specified proteome: please ensure that the specified proteome is the identical proteome used in the LCD-Composer analysis.')
        print('NOTE: Only UniProt proteomes are supported for GO-term analysis.')
        exit()
    elif len(valid_prots) < len(lcd_prots):
        print('\nWARNING: The following proteins could not be found in the specified proteome file.')
        print([prot for prot in lcd_prots if prot not in valid_prots])
        

def output_GOresults(random_id, fasta_str, g_res):

    output = open(random_id+'_GOterm_Results.tsv', 'w')
    new_header = ['# GO', 'NS', 'enrichment', 'name', 'hits_in_study', 'totalProteins_in_study', 'ratio_in_study', 'hits_in_population', 'totalProteins_in_population', 'ratio_in_pop', 'odds_ratio', 'p_uncorrected', 'depth', 'study_count', 'p_bonferroni', 'p_sidak', 'p_holm', 'study_items']
    if fasta_str in gaf_df:
        length_assoc_gos = get_length_assoc_gos(fasta_str)
        new_header.insert(-1, 'Length-associated GO term?')
    output.write('\t'.join(new_header) + '\n')

    lines = []
    sidak_pvals = []
    oddsratios = []
    for enr_record in g_res:
        items = str(enr_record).split('\t')
        go_id = items[0]

        ratio_in_study, ratio_in_pop = items[4:6]
        hits_study, total_study = ratio_in_study.split('/')
        hits_pop, total_pop = ratio_in_pop.split('/')
        sidak_pval = float(items[10])
        
        contingency = [[int(hits_study), (int(total_study)-int(hits_study))], [int(hits_pop), (int(total_pop)-int(hits_pop))]]
        
        oddsratio, pval = stats.fisher_exact(contingency, alternative='two-sided')
        
        new_line = [go_id] + items[1:4] + [hits_study, total_study] + [items[4]] + [hits_pop, total_pop] + [items[5]] + [str(oddsratio)] + items[6:]
        if fasta_str in gaf_df:
            if go_id in length_assoc_gos:
                new_line.insert(-1, '1')
            else:
                new_line.insert(-1, '0')
        lines.append(new_line)
        sidak_pvals.append(sidak_pval)
        oddsratios.append(oddsratio)
        
    # PERFORM SEQUENTIAL SORTING. SORTS BY oddsratio (HIGH-TO-LOW), THEN SORTS BY sidak_pval (LOW-TO-HIGH). END RESULT IS THAT LOWEST P-VALUES SHOW UP FIRST, THEN FOR ALL LINES WITH THE SAME P-VALUE, THOSE WITH THE HIGHEST oddsratio SHOW UP FIRST.
    oddsratios, lines, sidak_pvals = zip(*sorted(zip(oddsratios, lines, sidak_pvals), reverse=True))
    heirarchy = [x for x in range(len(sidak_pvals))]
    sidak_pvals, heirarchy, oddsratios, lines = zip(*sorted(zip(sidak_pvals, heirarchy, oddsratios, lines)))
    for line in lines:
        output.write('\t'.join(line) + '\n')

    output.close()
    
    
def get_length_assoc_gos(fasta_str):
        
    h = open('Length-Associated_GOterms_SupportedOrganisms.tsv')
    header = h.readline()
    
    length_assoc_gos = []
    for line in h:
        go_id, goterm, organism, proteome = line.rstrip().split('\t')
        if fasta_str == proteome:
            length_assoc_gos.append(go_id)
            
    h.close()
    
    return length_assoc_gos
        
def get_args(arguments):
    parser = argparse.ArgumentParser(description='Test for LCD enrichment among a pre-defined set of proteins.')
    
    parser.add_argument('jobfile_id', help="""Your sequence file (in FASTA format).""")
    
    parser.add_argument('fasta_file', help="""The FASTA file containing your proteome of interest. This will serve as the "background" set of proteins for calculating GO term enrichment/purification among your LCD-containing proteins.""")
    
    parser.add_argument('gaf_file', help="""The name of the appropriate gene annotation file (GAF). This should correspond to the proteome that you are using to define the "background" set of proteins.""")

    args = parser.parse_args(arguments)
    
    return args
    

if __name__ == '__main__':
    import sys, argparse, math
    from goatools import obo_parser
    from goatools.go_enrichment import GOEnrichmentStudy
    from scipy import stats
    args = get_args(sys.argv[1:])
    main(args)