
"""
Description:
    LCD-Composer is a composition-based method for identifying low-complexity domains in protein sequences.
    Refer to Cascarina et. al. (2021) (https://academic.oup.com/nargab/article/3/2/lqab048/6285187) for a complete 
    description of the original LCD-Composer algorithm, which is used in original and modified form here.
    
    This script is a modification/augmentation of the original LCD-Composer algorithm, developed as 
    part of the LCD-Composer webserver.
    
    For a description of modifications and added functionality, refer to Cascarina and Ross (2022) 
    (https://pubmed.ncbi.nlm.nih.gov/36282522/) and the LCD-Composer webserver HELP page at 
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


residues = 'ACDEFGHIKLMNPQRSTVWY'

def LCDcomposer(random_id, fasta_files, win_size, amino_acids, other_aas, comp_thresholds, disp_threshold, disp_method, ign_disp):

    all_aas = 'ACDEFGHIKLMNPQRSTVWY'
    combined_aas = ''.join(amino_acids)

    mins_df = {}
    maxs_df = {}
    
    print('Start time:', str(datetime.datetime.now()))
    
    job_file = str(random_id) + '_LCDcomposer_RESULTS_WholeProteome.tsv'
    output = open(job_file, 'w')
    output_params_header(output, random_id, fasta_files, win_size, amino_acids, other_aas, comp_thresholds, disp_threshold, disp_method, ign_disp)
    output.write('\t'.join( ['Protein Description','UniProt ID (when applicable)','Domain Sequence','Domain Boundaries','Final Domain Composition','Final Domain Linear Dispersion','\t'.join([aa for aa in all_aas]) ]) + '\n')

    lcd_prots = set()
    
    proteome_size = get_proteome_size(fasta_files[0])
    if len(fasta_files) > 1:
        proteome_size += get_proteome_size(fasta_files[1])
    
    total_seqs = 0
    perc_range = [int(proteome_size*x/10) for x in range(1, 10)]
    
    for file_index, file in enumerate(fasta_files):
        h = open(file)

        for id, seq in fasta_parser(h):

            #PRINT PROGRESS AT 10% MILESTONES
            if total_seqs in perc_range:
                print(str(round(total_seqs/proteome_size*100)) + '% Complete' + '\t' + str(datetime.datetime.now()))
            
            if id.count('|') == 2:
                junk, uniprot, junk = id.split('|')
            else:
                uniprot = 'N/A'
                
            #REMOVE STOP CODON FROM C-TERMINUS AND WARN USERS IF A SEQUENCE CONTAINS MULTIPLE STOP CODONS
            if seq.count('*') == 1 and seq[-1] == '*':
                seq = seq[:-1]
            elif seq.count('*') > 1:
                print('\nRuntime warning: protID \"' + id[:15] + '...\" contains multiple stop codons, which will slightly affect the amino acid composition and separation calculations. Stop codons are automatically removed by the program from the C-terminus of each sequence but internal stop codons are not removed. Consider removing extra stop codons before evaluating, removing these sequences from analyses entirely, or evaluating sequences as-is.\n')
                if seq[-1] == '*':
                    seq = seq[:-1]

            normed_stddevs = []
            group_comps = {aa_group:[] for aa_group in amino_acids}
            hit_positions = []
            pos = -1
            for i in range(len(seq) - win_size+1):
                window = seq[i:i+win_size]
                    
                pos += 1

                #CALCULATE COMPOSITION OF EACH aa_group FOR THE WINDOW
                comps, group_comps = calc_composition(window, amino_acids, group_comps)
                
                #MASK REGIONS THAT DON'T PASS COMPOSITION THRESHOLD
                if 0 in [1 if comps[i] > comp_thresholds[i] else 0 for i in range(len(amino_acids))]:
                    normed_stddevs.append( -1 )
                    continue

                #PERFORM FULL CALCULATIONS
                else:
                    #QUICK CALCULATION FOR WINDOWS THAT ARE 100% AAs OF INTEREST
                    if sum([window.count(aa)/win_size*100 for aa in combined_aas]) > 100-0.000001:
                        normed_stddevs.append(1)
                        hit_positions.append(pos)
                        continue
                        
                    norm_stddev = calc_dispersion(window, combined_aas, other_aas, mins_df, maxs_df)
                    normed_stddevs.append( norm_stddev )

                    if sum(comps) >= ign_disp or norm_stddev > disp_threshold:
                        hit_positions.append( pos )

            total_seqs += 1
            
            #DERIVE DOMAIN BOUNDARIES AND CORRESPONDING DOMAIN SEQUENCES
            domain_boundaries = merge_windows(hit_positions, win_size, id, seq)
            seqs = [seq[ bounds[0] : bounds[1] ] for bounds in domain_boundaries]

            #TRIM TERMINI UNTIL THEY MATCH THE RESIDUE OF INTEREST (OR ONE OF THE RESIDUES OF INTEREST FOR AA GROUPS)
            trimmed_seqs, trimmed_boundaries = trim_termini(seqs, combined_aas, domain_boundaries)

            #CONVERT TO MORE USER-INTUITIVE PROTEIN BOUNDARY NUMBERING
            trimmed_boundaries = ['('+str(bounds[0]+1) + '-' + str(bounds[1])+')' for bounds in trimmed_boundaries]
            final_comps, final_stddevs = calc_final_comps_stddevs(trimmed_seqs, combined_aas, win_size, mins_df, maxs_df, other_aas)
            
            junk, uniprot, *junk = id.split('|')
            if len(trimmed_boundaries) > 0:
                lcd_prots.add(uniprot)
                for i in range(len(trimmed_seqs)):
                    output.write('\t'.join( [id, uniprot, trimmed_seqs[i], trimmed_boundaries[i], str(round(final_comps[i], 2)), str(round(final_stddevs[i], 4)), '\t'.join([str(round(trimmed_seqs[i].count(aa) / len(trimmed_seqs[i]) * 100, 2)) for aa in all_aas]) ]) + '\n')
            else:
                continue
            
        h.close()
    output.close()
    
    print('LCD search complete:', str(datetime.datetime.now()))
    print('Starting LCD enrichment analysis...')
    
    return job_file, lcd_prots
    
    
def LCDcomposer_NewDispMethod(random_id, fasta_files, win_size, amino_acids, other_aas, comp_thresholds, disp_thresholds, disp_method, ign_disp):

    all_aas = 'ACDEFGHIKLMNPQRSTVWY'
    combined_aas = ''.join(amino_acids)
    
    other_aas_groupform = []
    for aa_group in amino_acids:
        other_groupform = all_aas[:]
        for aa in aa_group:
            other_groupform = other_groupform.replace(aa, '')
        other_aas_groupform.append(other_groupform)

    mins_df = {}
    maxs_df = {}
    
    print('Start time:', str(datetime.datetime.now()))
    
    job_file = str(random_id) + '_LCDcomposer_RESULTS_WholeProteome.tsv'
    output = open(job_file, 'w')
    output_params_header(output, random_id, fasta_files, win_size, amino_acids, other_aas, comp_thresholds, disp_thresholds, disp_method, ign_disp)
    output.write('\t'.join( ['Protein Description','UniProt ID (when applicable)','Domain Sequence','Domain Boundaries','Final Domain Composition','Final Domain Linear Dispersion','\t'.join([aa for aa in all_aas]) ]) + '\n')

    lcd_prots = set()
    
    proteome_size = get_proteome_size(fasta_files[0])
    if len(fasta_files) > 1:
        proteome_size += get_proteome_size(fasta_files[1])
    
    total_seqs = 0
    perc_range = [int(proteome_size*x/10) for x in range(1, 10)]
    
    for file_index, file in enumerate(fasta_files):
        h = open(file)

        for id, seq in fasta_parser(h):

            #PRINT PROGRESS AT 10% MILESTONES
            if total_seqs in perc_range:
                print(str(round(total_seqs/proteome_size*100)) + '% Complete' + '\t' + str(datetime.datetime.now()))
            
            if id.count('|') == 2:
                junk, uniprot, junk = id.split('|')
            else:
                uniprot = 'N/A'
                
            #REMOVE STOP CODON FROM C-TERMINUS AND WARN USERS IF A SEQUENCE CONTAINS MULTIPLE STOP CODONS
            if seq.count('*') == 1 and seq[-1] == '*':
                seq = seq[:-1]
            elif seq.count('*') > 1:
                print('\nRuntime warning: protID \"' + id[:15] + '...\" contains multiple stop codons, which will slightly affect the amino acid composition and separation calculations. Stop codons are automatically removed by the program from the C-terminus of each sequence but internal stop codons are not removed. Consider removing extra stop codons before evaluating, removing these sequences from analyses entirely, or evaluating sequences as-is.\n')
                if seq[-1] == '*':
                    seq = seq[:-1]

            normed_stddevs = []
            group_comps = {aa_group:[] for aa_group in amino_acids}
            hit_positions = []
            pos = -1
            for i in range(len(seq) - win_size+1):
                window = seq[i:i+win_size]
                    
                pos += 1

                #CALCULATE COMPOSITION OF EACH aa_group FOR THE WINDOW
                comps, group_comps = calc_composition(window, amino_acids, group_comps)
                
                #MASK REGIONS THAT DON'T PASS COMPOSITION THRESHOLD
                if 0 in [1 if comps[i] > comp_thresholds[i] else 0 for i in range(len(amino_acids))]:
                    normed_stddevs.append( -1 )
                    continue
                elif sum(comps) < ign_disp and 0 in [1 if calc_dispersion(window, amino_acids[i], other_aas_groupform[i], mins_df, maxs_df) > disp_thresholds[i] else 0 for i in range(len(amino_acids))]:
                    normed_stddevs.append( -1 )
                    continue

                #PERFORM FULL CALCULATIONS
                else:
                    #QUICK CALCULATION FOR WINDOWS THAT ARE 100% AAs OF INTEREST
                    if sum([window.count(aa)/win_size*100 for aa in combined_aas]) > 100-0.000001:
                        normed_stddevs.append(1)
                        hit_positions.append(pos)
                        continue
                        
                    norm_stddev_vals = [calc_dispersion(window, amino_acids[i], other_aas_groupform[i], mins_df, maxs_df) for i in range(len(amino_acids))]
                    normed_stddevs.append( '_'.join([str(x) for x in norm_stddev_vals]) )
                    disp_checksum = sum( [1 if norm_stddev_vals[i] > disp_thresholds[i] else 0 for i in range(len(disp_thresholds))] )
                    
                    if sum(comps) >= ign_disp or disp_checksum == len(disp_thresholds):
                        hit_positions.append( pos )

            total_seqs += 1
            
            #DERIVE DOMAIN BOUNDARIES AND CORRESPONDING DOMAIN SEQUENCES
            domain_boundaries = merge_windows(hit_positions, win_size, id, seq)
            seqs = [seq[ bounds[0] : bounds[1] ] for bounds in domain_boundaries]

            #TRIM TERMINI UNTIL THEY MATCH THE RESIDUE OF INTEREST (OR ONE OF THE RESIDUES OF INTEREST FOR AA GROUPS)
            trimmed_seqs, trimmed_boundaries = trim_termini(seqs, combined_aas, domain_boundaries)

            #CONVERT TO MORE USER-INTUITIVE PROTEIN BOUNDARY NUMBERING
            trimmed_boundaries = ['('+str(bounds[0]+1) + '-' + str(bounds[1])+')' for bounds in trimmed_boundaries]
            final_comps, final_stddevs = calc_final_comps_stddevs_NewMethod(trimmed_seqs, combined_aas, win_size, mins_df, maxs_df, other_aas, amino_acids, other_aas_groupform)
            
            junk, uniprot, *junk = id.split('|')
            if len(trimmed_boundaries) > 0:
                lcd_prots.add(uniprot)
                for i in range(len(trimmed_seqs)):
                    output.write('\t'.join( [id, uniprot, trimmed_seqs[i], trimmed_boundaries[i], str(round(final_comps[i], 2)), final_stddevs[i], '\t'.join([str(round(trimmed_seqs[i].count(aa) / len(trimmed_seqs[i]) * 100, 2)) for aa in all_aas]) ]) + '\n')
            else:
                continue
            
        h.close()
    output.close()
    
    print('LCD search complete:', str(datetime.datetime.now()))
    print('Starting LCD enrichment analysis...')
    
    return job_file, lcd_prots
    
    
def LCDcomposer_AutoDetect(random_id, fasta_files, win_size, amino_acids, other_aas, comp_thresholds, disp_threshold, disp_method, ign_disp):

    all_aas = 'ACDEFGHIKLMNPQRSTVWY'

    mins_df = {}
    maxs_df = {}
    
    print('Start time:', str(datetime.datetime.now()))
    
    job_file = str(random_id) + '_LCDcomposer_RESULTS_WholeProteome.tsv'
    output = open(job_file, 'w')
    output_params_header(output, random_id, fasta_files, win_size, amino_acids, other_aas, comp_thresholds, disp_threshold, disp_method, ign_disp)
    output.write('\t'.join( ['Protein Description','UniProt ID (when applicable)','Domain Sequence','Domain Boundaries','LCD Class','Final Domain Composition','Final Domain Linear Dispersion','\t'.join([aa for aa in all_aas]) ]) + '\n')

    lcd_prots = {aa:set() for aa in all_aas}
    
    proteome_size = get_proteome_size(fasta_files[0])
    if len(fasta_files) > 1:
        proteome_size += get_proteome_size(fasta_files[1])
    
    total_seqs = 0
    perc_range = [int(proteome_size*x/10) for x in range(1, 10)]
    
    for file_index, file in enumerate(fasta_files):
        h = open(file)

        for id, seq in fasta_parser(h):

            #PRINT PROGRESS AT 10% MILESTONES
            if total_seqs in perc_range:
                print(str(round(total_seqs/proteome_size*100)) + '% Complete' + '\t' + str(datetime.datetime.now()))
            
            if id.count('|') == 2:
                junk, uniprot, junk = id.split('|')
            else:
                uniprot = 'N/A'
                
            #REMOVE STOP CODON FROM C-TERMINUS AND WARN USERS IF A SEQUENCE CONTAINS MULTIPLE STOP CODONS
            if seq.count('*') == 1 and seq[-1] == '*':
                seq = seq[:-1]
            elif seq.count('*') > 1:
                print('\nRuntime warning: protID \"' + id[:15] + '...\" contains multiple stop codons, which will slightly affect the amino acid composition and separation calculations. Stop codons are automatically removed by the program from the C-terminus of each sequence but internal stop codons are not removed. Consider removing extra stop codons before evaluating, removing these sequences from analyses entirely, or evaluating sequences as-is.\n')
                if seq[-1] == '*':
                    seq = seq[:-1]
            
            for res in 'ACDEFGHIKLMNPQRSTVWY':
                amino_acids = [res]     # OVERWRITES amino_acids EACH TIME. THIS IS SPECIFIC FOR THE AUTO-DETECT MODE.
                combined_aas = res[:]
                other_aas = all_aas.replace(res, '')
                normed_stddevs = []
                group_comps = {aa_group:[] for aa_group in res}
                hit_positions = []
                pos = -1
                for i in range(len(seq) - win_size+1):
                    window = seq[i:i+win_size]
                        
                    pos += 1

                    #CALCULATE COMPOSITION OF EACH aa_group FOR THE WINDOW
                    comps, group_comps = calc_composition(window, amino_acids, group_comps)
                    
                    #MASK REGIONS THAT DON'T PASS COMPOSITION THRESHOLD
                    if 0 in [1 if comps[i] > comp_thresholds[i] else 0 for i in range(len(amino_acids))]:
                        normed_stddevs.append( -1 )
                        continue

                    #PERFORM FULL CALCULATIONS
                    else:
                        #QUICK CALCULATION FOR WINDOWS THAT ARE 100% AAs OF INTEREST
                        if sum([window.count(aa)/win_size*100 for aa in combined_aas]) > 100-0.000001:
                            normed_stddevs.append(1)
                            hit_positions.append(pos)
                            continue
                            
                        norm_stddev = calc_dispersion(window, combined_aas, other_aas, mins_df, maxs_df)
                        normed_stddevs.append( norm_stddev )
                        if sum(comps) >= ign_disp or norm_stddev > disp_threshold:
                            hit_positions.append( pos )
                
                #DERIVE DOMAIN BOUNDARIES AND CORRESPONDING DOMAIN SEQUENCES
                domain_boundaries = merge_windows(hit_positions, win_size, id, seq)
                seqs = [seq[ bounds[0] : bounds[1] ] for bounds in domain_boundaries]

                #TRIM TERMINI UNTIL THEY MATCH THE RESIDUE OF INTEREST (OR ONE OF THE RESIDUES OF INTEREST FOR AA GROUPS)
                trimmed_seqs, trimmed_boundaries = trim_termini(seqs, combined_aas, domain_boundaries)

                #CONVERT TO MORE USER-INTUITIVE PROTEIN BOUNDARY NUMBERING
                trimmed_boundaries = ['('+str(bounds[0]+1) + '-' + str(bounds[1])+')' for bounds in trimmed_boundaries]
                final_comps, final_stddevs = calc_final_comps_stddevs(trimmed_seqs, combined_aas, win_size, mins_df, maxs_df, other_aas)
                
                junk, uniprot, *junk = id.split('|')
                if len(trimmed_boundaries) > 0:
                    lcd_prots[res].add(uniprot)
                    for i in range(len(trimmed_seqs)):
                        output.write('\t'.join( [id, uniprot, trimmed_seqs[i], trimmed_boundaries[i], res, str(round(final_comps[i], 2)), str(round(final_stddevs[i], 4)), '\t'.join([str(round(trimmed_seqs[i].count(aa) / len(trimmed_seqs[i]) * 100, 2)) for aa in all_aas]) ]) + '\n')
                else:
                    continue
        
            total_seqs += 1
            
        h.close()

    output.close()
    
    print('LCD search complete:', str(datetime.datetime.now()))
    print('Starting LCD enrichment analysis...')
    
    return job_file, lcd_prots
    

def calc_composition(window, amino_acids, group_comps):
    """
    Calculates the amino acid composition for each user-defined amino acid group in amino_acids.
    
    Returns:
        Compositions (list)
        Group compositions (dictionary)
    """
    comps = []
    for aa_group in amino_acids:
        comp = 0
        for aa in aa_group:
            comp += window.count(aa) / len(window) * 100
        comps.append( comp )
        group_comps[aa_group].append( comp )
        
    return comps, group_comps
    

def calc_dispersion(window, combined_aas, other_aas, mins_df, maxs_df):
    """
    Calculates the normalized spacing standard deviation for a given window sequence.
    
    Returns:
        Normalized spacing standard deviations (float)
    """
    
    res_spacing_stddev = calc_spacing(window, combined_aas, other_aas)
    num_muts = sum([window.count(aa) for aa in other_aas])

    min_stddev, max_stddev, mins_df, maxs_df = get_MinMax_stddevs(len(window), num_muts, mins_df, maxs_df)
    
    #HANDLES CASES WITH 100% COMPOSITION (min_stddev AND max_stddev ARE BOTH ZERO)
    if min_stddev == max_stddev:
        norm_stddev = 1.0
    else:
        norm_stddev = 1 - (res_spacing_stddev - min_stddev) / (max_stddev - min_stddev)
        
    return norm_stddev
                    

def calc_final_comps_stddevs(trimmed_seqs, combined_aas, win_size, mins_df, maxs_df, other_aas):
    """
    Calculates the final composition and normalized spacing standard deviation 
    for each merged & trimmed domain identified.
    
    Returns:
        Compositions (list)
        Normalized spacing standard deviations (list)
    """
    
    final_comps = []
    final_stddevs = []
    for seq in trimmed_seqs:
        comp = sum([seq.count(aa) for aa in combined_aas]) / len(seq) * 100
        norm_stddev = calc_dispersion(seq, combined_aas, other_aas, mins_df, maxs_df)
        num_muts = sum([seq.count(aa) for aa in other_aas])
        
        final_comps.append(comp)
        final_stddevs.append(norm_stddev)
        
    return final_comps, final_stddevs
    
def calc_final_comps_stddevs_NewMethod(trimmed_seqs, combined_aas, win_size, mins_df, maxs_df, other_aas, amino_acids, other_aas_groupform):
    """
    Calculates the final composition and normalized spacing standard deviation 
    for each merged & trimmed domain identified.
    
    Returns:
        Compositions (list)
        Normalized spacing standard deviations (list)
    """
    
    final_comps = []
    final_stddevs = []
    for seq in trimmed_seqs:
        comp = sum([seq.count(aa) for aa in combined_aas]) / len(seq) * 100
        norm_stddevs = [calc_dispersion(seq, amino_acids[i], other_aas_groupform[i], mins_df, maxs_df) for i in range(len(amino_acids)) ]
        num_muts = sum([seq.count(aa) for aa in other_aas])
        
        final_comps.append(comp)
        final_stddevs.append('_'.join([str(round(x, 4)) for x in norm_stddevs]))
        
    return final_comps, final_stddevs
    

def get_MinMax_stddevs(win_size, num_muts, mins_df, maxs_df):
    """
    Checks whether the minimum possible stddev and maximum possible stddev have been 
    calculated for a sequence of length win_size and a composition 
    corresponding to num_muts.
    
    If they have not already been calculated, generates sequences with minimum and 
    maximum standard deviations in spacings and adds them to lookup dictionaries.
    
    Returns:
        min_stddev (float)
        max_stddev (float)
        mins_df (dictionary)
        maxs_df (dictionary)
    """

    mins_df[win_size] = mins_df.get(win_size, {})
    maxs_df[win_size] = maxs_df.get(win_size, {})
    
    if num_muts not in mins_df[win_size]:
        min_seq = generate_min_seq(win_size, num_muts)
        max_seq = generate_max_seq(win_size, num_muts)

        min_stddev = calc_spacing(min_seq, 'Z', 'G')
        max_stddev = calc_spacing(max_seq, 'Z', 'G')

        mins_df[win_size][num_muts] = mins_df[win_size].get(num_muts, min_stddev)
        maxs_df[win_size][num_muts] = maxs_df[win_size].get(num_muts, max_stddev)

    return mins_df[win_size][num_muts], maxs_df[win_size][num_muts], mins_df, maxs_df
                    
                    
def trim_termini(seqs, amino_acids, domain_boundaries):
    """
    Trims domains that pass composition/separation thresholds until the termini match 
    one of the specified amino acids of interest.
    
    Returns:
        Trimmed sequences (list).
        Trimmed domains (list of tuples, where each tuple contains the start and 
                        end indices of the domain, respectively).
    """
    
    for i in range(len(seqs)):
        while seqs[i][0] not in amino_acids:
            seqs[i] = seqs[i][1:]
            domain_boundaries[i] = (domain_boundaries[i][0]+1, domain_boundaries[i][1])
            
        while seqs[i][-1] not in amino_acids:
            seqs[i] = seqs[i][:-1]
            domain_boundaries[i] = (domain_boundaries[i][0], domain_boundaries[i][1]-1)

    return seqs, domain_boundaries
    

def merge_windows(hit_positions, win_size, id, seq):
    """
    Merge neighboring windows that pass the composition/separation thresholds and
    lie within 1/2 the window size distance from each other.
    
    Returns:
        The boundaries (according to Python indexing, not protein positions) of all identified domains.
        This is a list of tuples, where each tuple contains the start and 
        end indices of the domain, respectively.
    """

    j = 0
    domain_boundaries = []
    
    while hit_positions:
        start = hit_positions[0]
        while j < (len(hit_positions) - 1) and ( (hit_positions[j] + win_size) >= ( hit_positions[j+1] ) ):
            j += 1
        
        #GETS THE ENDING POSITION FOR THE CURRENT WINDOW AND STORES THE 2-TUPLE
        end = hit_positions[j]
        domain_boundaries.append( (start , end+win_size) )

        #MODIFIES hit_positions TO DELETE ALL POSITIONS THAT WERE JUST MERGED, THEN RE-SETS j=0 TO START ITERATING AT THE FIRST POSITION IN THE NEW LIST
        hit_positions = hit_positions[j+1:]
        j = 0

    return domain_boundaries
    
        
def calc_spacing(seq, amino_acids, other_aas):
    """
    Calculates the spacings of the either the amino acid(s) of interest or all other amino acids (whichever is smallest).
    
    Returns:
        Standard deviation of the spacings (float)
    """

    #ADDS RESIDUE OF INTEREST TO THE ENDS SO THAT SPACINGS ARE ALSO CALCULATED BETWEEN THE N-TERM AND THE FIRST RESIDUE OF INTEREST, AND BETWEEN THE LAST RESIDUE OF INTEREST AND THE C-TERM
    amino_acids = ''.join(amino_acids)
    seq = amino_acids[0] + seq + amino_acids[0]
    res_pos = [x for x in range(len(seq)) if seq[x] in amino_acids]
    
    seq = other_aas[0] + seq[1:-1] + other_aas[0]
    nonres_pos = [x for x in range(len(seq)) if seq[x] not in amino_acids]
    res_spacings = sorted( [ res_pos[x] - res_pos[x-1] for x in range(1, len(res_pos))] +  [ nonres_pos[x] - nonres_pos[x-1] for x in range(1, len(nonres_pos))])
    
    return calc_stddev(res_spacings)
    
    
def calc_stddev(spacings):
    """
    Calculates the population standard deviation.
    
    Returns:
        Standard deviation (float)
    """
    
    mean = statistics.mean(spacings)
    var = sum( [(x - mean)**2 for x in spacings] ) / len(spacings)
    stddev = math.sqrt(var)
    
    return stddev
    
        
def get_proteome_size(fasta_file):
    """Function to quickly tally the # of proteins in the user-defined FASTA file
    
    Returns:
        Proteome size (int)
    """
    h = open(fasta_file)
    proteome_size = 0
    for line in h:
        if line.startswith('>'):
            proteome_size += 1
            
    return proteome_size
    
    
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
                

def generate_min_seq(win_size, num_muts):
    """
    Generate a sequence with Z num_muts and maximum spacing between Z's.
    This will be used to calculate the minimum stddev value for normalization,
    since large spacing results in small stddevs.
    
    Returns:
        min_seq (string)
    """
    
    #CORNER CASE WITH SINGLE MUTATION
    if num_muts == 1:
        min_seq = 'G'*int(win_size/2) + 'Z'
        while len(min_seq) < win_size:
            min_seq += 'G'
            
    #CORNER CASES WHERE NUMBER OF MUTATIONS IS ~1/2 THE win_size
    elif num_muts == int(win_size / 2):
        num_muts = int(win_size / 2) + 1
        min_seq = 'ZG'*num_muts
            
    else:
        if num_muts >= int(win_size/2):
            num_muts = win_size-num_muts
            
        min_seq = ''
        block_size = int( (win_size - num_muts) / (num_muts+1) )
        remainder = (win_size-num_muts) % (num_muts+1)
        
        while remainder > 0:
            min_seq += 'Z'*(block_size+1) + 'G'
            remainder -= 1
            
        while len(min_seq) < win_size:
            min_seq += 'Z'*block_size + 'G'
            
    min_seq = min_seq[:win_size]
    
    return min_seq
    
    
def generate_max_seq(win_size, num_muts):
    """
    Generate a sequence with Z num_muts and all Z's clustered at one end (minimum spacing)
    This will be used to calculate the maximum stddev value for normalization.
    
    Returns:
        max_seq (string)
    """
    
    max_seq = num_muts*'Z' + 'G'*(win_size-num_muts)
    
    return max_seq
    
    
def holmsidak_correction(pvals):
    
    positions = [i for i in range(len(pvals))]
    sorted_pvals, pval_positions = zip(*sorted(zip(pvals, positions)))
    m = len(pvals)

    mpmath.mp.dps = 300

    corrected_pvals = [1-(1-mpmath.mp.mpf(sorted_pvals[len(pvals)-i]))**(i) for i in range(m, 0, -1)]
    
    final_pvals = [corrected_pvals[0]]
    for i in range(1, len(corrected_pvals)):
        pval = max(final_pvals[-1], corrected_pvals[i])
        final_pvals.append(pval)

    final_pvals = [mpmath.nstr(x, 16) for x in final_pvals]

    pval_positions, final_pvals = zip(*sorted(zip(pval_positions, final_pvals)))

    return final_pvals
    

def crosscheck_user_prots(user_uniprots, fasta_files):
    
    proteome_prots = set()
    no_isoforms_prots = set()
    for file_index, file in enumerate(fasta_files):
        h = open(file)
        if file_index == 0:     #MILDLY REPETITIVE CODE STRUCTURE BUT SHOULD BE SLIGHTLY FASTER THAN RUNNING "if" CHECK WITHIN THE FOR LOOP.
            for id, seq in fasta_parser(h):
                junk, uniprot, *junk = id.split('|')
                proteome_prots.add(uniprot)
                no_isoforms_prots.add(uniprot)
        else:
            for id, seq in fasta_parser(h):
                junk, uniprot, *junk = id.split('|')
                proteome_prots.add(uniprot)
            
        h.close()
        
    missing_prots = set()
    all_prots_for_GOanalysis = set()
    for prot in user_uniprots:
        if prot not in proteome_prots:
            missing_prots.add(prot)
        if prot in no_isoforms_prots:
            all_prots_for_GOanalysis.add(prot)
            
    if len(missing_prots) == len(user_uniprots):
        print('\nError:\nNone of the UniProt identifiers could be found in the organism selected.\n')
        exit()
            
    return proteome_prots, missing_prots, all_prots_for_GOanalysis, no_isoforms_prots
    
    
def get_LCDprots(job_file):
    
    h = open(job_file)
    for i in range(11):
        h.readline()
        
    lcd_prots = set()
    for line in h:
        fasta_header, uniprot, *junk = line.rstrip().split('\t')
        lcd_prots.add(uniprot)
    h.close()
    
    return lcd_prots
    
    
def get_proteome_prots(fasta_files, user_uniprots):

    prots = set()
    for file in fasta_files:
        h = open(file)
        for id, seq in fasta_parser(h):
            junk, uniprot, *junk = id.split('|')
            if uniprot not in user_uniprots:
                prots.add(uniprot)
            
        h.close()
    
    return prots
    
    
def test_LCD_enrichment(lcd_prots, user_uniprots, proteome_prots, missing_prots, filtered_proteome_prots, impute_value):
    
    user_uniprots = {prot for prot in user_uniprots if prot not in missing_prots}
    user_hits = 0
    user_lcdprots = set()
    for user_prot in user_uniprots:
        if user_prot in lcd_prots:
            user_hits += 1
            user_lcdprots.add(user_prot)
            
    proteome_hits = 0
    proteome_lcdprots = set()
    for prot in lcd_prots:
        if prot in filtered_proteome_prots:
            proteome_hits += 1
            proteome_lcdprots.add(prot)
            
    user_nonhits = len(user_uniprots) - user_hits
    proteome_nonhits = len(filtered_proteome_prots) - proteome_hits
    imputation_used = False
    contingency = [[user_hits, user_nonhits], [proteome_hits, proteome_nonhits]]
    if user_hits == 0 and impute_value and proteome_hits != 0:
        oddsratio, pval = stats.fisher_exact(contingency, alternative='two-sided')
        contingency = [[1, user_nonhits], [proteome_hits, proteome_nonhits]]
        oddsratio, biased_pval = stats.fisher_exact(contingency, alternative='two-sided')
        imputation_used = True
    elif proteome_hits == 0 and impute_value and user_hits != 0:
        oddsratio, pval = stats.fisher_exact(contingency, alternative='two-sided')
        contingency = [[user_hits, user_nonhits], [1, proteome_nonhits]]
        oddsratio, biased_pval = stats.fisher_exact(contingency, alternative='two-sided')
        imputation_used = True
    else:
        oddsratio, pval = stats.fisher_exact(contingency, alternative='two-sided')

    try:
        log2OR = math.log2(oddsratio)
    except:
        log2OR = -1000    # DUMMY VALUE SO THAT SUBSEQUENT SORTING OPERATION STILL WORKS. THIS IS LATER RELACED WITH "N/A" STRING.
    
    return user_hits, user_nonhits, proteome_hits, proteome_nonhits, log2OR, pval, user_lcdprots, proteome_lcdprots, imputation_used
    
    
def write_statsfile(random_id, rows, header_labels):

    filename  = random_id + '_LCD-EnrichmentTest_Statistics.tsv'
    output = open(filename, 'w')
    header_labels += ['IDs of Your Proteins with an LCD', 'IDs of Proteome Proteins with an LCD (excluding your proteins)']
    output.write('\t'.join(header_labels) + '\n')
    for row in rows:
        output.write('\t'.join( [str(x) for x in row] ) + '\n')
        
    output.close()
    
    return filename
    

def barplot(random_id, df, color_palette):

    sns.barplot(x='Amino Acid', y='log2OR', data=df, palette=color_palette)
    pvals = df['P-value']

    ax = plt.gca()
    rects = [rect for rect in ax.get_children() if isinstance(rect, patches.Rectangle)]

    offset = 0.04
    heights = []
    for i in range(len(rects)-1):
        if i > len(pvals):
            continue
        height = rects[i].get_height()
        heights.append(height)

        pval = pvals[i]

        if pval < 0.001:
            pval_indicator = '***'
        elif pval < 0.01:
            pval_indicator = '**'
        elif pval < 0.05:
            pval_indicator = '*'
        else:
            pval_indicator = ''

        if pval_indicator:
            if math.isinf(height):
                pval_indicator = '(inf)' + pval_indicator
            if height < 0:
                ax.annotate(pval_indicator, xy=(rects[i].get_x() + rects[i].get_width()/2., height), xytext=(0, -2.5), textcoords='offset points', ha='center', va='top')
            else:
                ax.annotate(pval_indicator, xy=(rects[i].get_x() + rects[i].get_width()/2., height), xytext=(0, 0), textcoords='offset points', ha='center')
        elif math.isinf(height):
            ax.annotate('(inf)', xy=(0,0.25), xytext=(0,0), textcoords='offset points', ha='center', rotation=90)

    plt.ylim(min(min([x for x in heights if not math.isinf(x)])-0.3, 0), max([x for x in heights if not math.isinf(x)])+0.5)
    plt.xticks(fontname='Arial', fontsize=14)
    plt.yticks(fontname='Arial', fontsize=14)
    plt.xlabel('LCD Category', fontname='Arial', fontsize=16)
    plt.ylabel('Degree of Enrichment', fontname='Arial', fontsize=16)
    fig = plt.gcf()
    fig.set_size_inches(8, 4.5)
    filename = random_id + '_LCD_Enrichment_Barplot.tiff'
    plt.savefig(filename, bbox_inches='tight', dpi=600, pil_kwargs={'compression':'tiff_lzw'})
    plt.close()
    
    return filename
    
    
def write_LCDresults_SubmittedProteins(random_id, results_file, fasta_files, win_size, amino_acids, other_aas, comp_thresholds, disp_threshold, disp_method, ign_disp, user_lcdprots):
    
    results_filename = str(random_id) + '_LCDcomposer_RESULTS.tsv'
    output = open(results_filename, 'w')
    output = output_params_header(output, random_id, fasta_files, win_size, amino_acids, other_aas, comp_thresholds, disp_threshold, disp_method, ign_disp)
    h = open(results_file)
    for i in range(10):
        h.readline()
        
    header = h.readline()
    output.write(header)
        
    for line in h:
        prot_desc, prot_id, *remainder = line.rstrip().split('\t')
        if prot_id not in user_lcdprots:
            continue
            
        output.write(line)
        
    h.close()
    output.close()
    
    return results_filename
    

def get_params(args):
    """Gather, define, and validate user parameters."""

    allowed_seq_chars = 'ACDEFGHIKLMNPQRSTVWY*'
    fasta_files = [args.fasta_file]
    isoforms_file = args.isoforms_file
    if isoforms_file:
        fasta_files.append(isoforms_file)
        
    allowed_uniprot_chars = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    allowed_uniprot_chars += '-1234567890'

    uniprots_file = args.uniprots_file
    try:
        h = open(uniprots_file)
        uniprots_list = []
        for line in h:
            uniprot = line.strip()
            uniprots_list.append(uniprot)
        h.close()

        if len(uniprots_list) > 0:
            char_check_str = ''.join( uniprots_list )
            if sum([1 for char in char_check_str if char in allowed_uniprot_chars]) != len(char_check_str):
                print('\nError:\nYour UniProt identifiers can only contain letters, numbers, and dash ("-") symbols.\n')
                exit()
        else:
            print('\nError:\nNo UniProt identifiers could be found in your provided uniprots file.\n')
            exit()
    except:
        print('\nError:\nThere was an error reading the file containing your UniProt identifiers. Please ensure that the file is correctly named and is in the same folder/directory as this script.\n')
        exit()
        
    run_autodetect = args.run_autodetect
    amino_acids = args.amino_acids
    if run_autodetect or amino_acids == 'Auto-detect':  
        amino_acids = 'Q'   # THIS IS A TEMPORARY VALUE WHEN auto-detect MODE IS USED TO ENSURE THAT ALL DOWNSTREAM PARAMETER CHECKS RUN PROPERLY.
        run_autodetect = True
    else:
        amino_acids = args.amino_acids.split('_')
        amino_acids = [x.upper() for x in amino_acids]
        aa_str = ''.join(amino_acids)
        if 0 in [0 for res in aa_str if res not in allowed_seq_chars[:-1] + '_']:
            print('\nError:\nThe "Amino Acid(s) of Interest" can only be composed of the 20 canonical amino acids, with underscores to separate amino acid groups.\n')
            exit()
    
    win_size = args.window_size
    try:
        win_size = int(win_size)
        if win_size < 1:
            exit()
    except:
        print('\nError:\nInvalid window size. The window size must be a positive integer.\n')
        exit()
        
    chars = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ*_'
    other_aas = ''.join([x for x in chars if x not in ''.join(amino_acids)])
    comp_thresholds = args.composition.split('_')
    comp_thresholds = [float(x)-0.000001 for x in comp_thresholds]    # +/-0.000001 --> ADJUST VALUE FOR DOWNSTREAM FLOATING POINT COMPARISONS
    try:
        comp_thresholds = [float(x)-0.000001 for x in comp_thresholds]    # +/-0.000001 --> ADJUST VALUE FOR DOWNSTREAM FLOATING POINT COMPARISONS
        if 0 in [1 if 0.999999<=comp_threshold<=100.0 else 0 for comp_threshold in comp_thresholds]:
            print('\nError:\nInvalid composition threshold(s). The composition threshold(s) must be values between 1-100 (inclusive).\n')
            exit()
    except:
        print('\nError:\nInvalid composition threshold(s). The composition threshold(s) must be values between 1-100 (inclusive).\n')
        exit()
        
    disp_threshold = args.dispersion
    
    disp_method = args.dispersion_method
    if disp_method == 'New' and not run_autodetect:
        if disp_threshold == 'default':
            disp_threshold = [0.5 for i in range(len(amino_acids))]
        else:
            try:
                disp_threshold = disp_threshold.split('_')
                disp_threshold = [float(x) - 0.0000001 for x in disp_threshold]    # +/-0.000001 --> ADJUST VALUE FOR DOWNSTREAM FLOATING POINT COMPARISONS
                threshold_check = [1 if (disp_threshold[i] > 0-0.000001 and disp_threshold[i] < 1+0.000001) else 0 for i in range(len(disp_threshold))]
                if 0 in threshold_check:
                    print('\nError:\nInvalid linear dispersion thresholds. Each linear dispersion threshold must be a number between 0.0-1.0 (inclusive)\n')
                    exit()
                if len(disp_threshold) != len(comp_thresholds):
                    print('\nError:\nWhen using the "New" method for linear dispersion, the number of linear dispersion thresholds must match the number of amino acid groups and composition thresholds you specify, each separated by an underscore (e.g. Amino acids: QN_G_Y, Composition thresholds: 40_15_10, Linear dispersion: 0.5_0.7_0.35).\n')
                    exit()
            except:
                print('\nError:\nThe linear dispersion threshold(s) must be numbers between 0.0-1.0 (inclusive).\n')
                exit()
    elif disp_method == 'Original' or run_autodetect:
        if disp_threshold == 'default':
            disp_threshold = 0.5
        else:
            try:
                disp_threshold = float(args.dispersion) - 0.000001    # +/-0.000001 --> ADJUST VALUE FOR DOWNSTREAM FLOATING POINT COMPARISONS
                if disp_threshold < 0-0.000001 or disp_threshold > 1+0.000001:
                    print('\nError:\nInvalid linear dispersion threshold. The linear dispersion threshold must be a number between 0.0-1.0 (inclusive)\n')
                    exit()
            except:
                print('\nError:\nThe linear dispersion threshold must be a number between 0.0-1.0 (inclusive).\n')
                exit()
    else:
        print('\nError:\nThe only acceptable arguments for the disp_method parameter are "New" and "Original".\n')
        exit()
    
    if not args.ignore_dispersion_threshold:
        ign_disp = sum(comp_thresholds) + ( (100 - sum(comp_thresholds)) / 2 ) - 0.000001
    else:
        try:
            ign_disp = float(args.ignore_dispersion_threshold) - 0.000001
        except:
            print('\nError:\nInvalid ignore dispersion threshold. Ignore dispersion threshold must be a decimal between 0.0-1.0.\n')
            exit()

    #RUN GATEKEEPER CHECKS=========================
    if sum(comp_thresholds) < 0-0.000001 or sum(comp_thresholds) > 100+0.000001:
        print('\nError:\nInvalid composition threshold. The composition threshold must be a number between 0-100 (inclusive)\n')
        exit()
        
    if len(amino_acids) != len(comp_thresholds) and not run_autodetect:
        print('\nError:\nInvalid composition thresholds for the specified amino acids or amino acid groups. The number of composition thresholds must match the number of amino acids or amino acid groups specified.\n')
        exit()

    if ign_disp < 0-0.000001 or ign_disp > 100+0.000001:
        print('\nError:\nInvalid ignore dispersion threshold. The ignore dispersion threshold must be a number between 0-100 (inclusive)\n')
        exit()
        
    if run_autodetect and len(comp_thresholds) > 1:
        print('\nError:\nInvalid composition threshold for the auto-detect mode. You may only specify a single composition threshold and a single linear dispersion threshold when using auto-detect mode.\n')
        exit()

    if run_autodetect:
        amino_acids = 'Auto-detect'

    impute_value = args.impute_value
    
    return fasta_files, win_size, amino_acids, other_aas, comp_thresholds, disp_threshold, disp_method, ign_disp, uniprots_list, run_autodetect, impute_value
    

def output_params_header(output, random_id, fasta_files, win_size, amino_acids, other_aas, comp_thresholds, disp_threshold, disp_method, ign_disp):

    output.write('>*RUNTIME PARAMETERS*\n')
    output.write('>Job ID: ' + str(random_id) + '\n')
    output.write('>FASTA File(s): ' + ', '.join(fasta_files) + '\n')
    output.write('>Window Size: ' + str(win_size) + '\n')
    output.write('>Amino Acid(s): ' + str(amino_acids) + '\n')
    output.write('>Composition Threshold(s): ' + str([round(x, 1) for x in comp_thresholds]) + '\n')
    try:
        output.write('>Linear Dispersion Threshold: ' + str([round(x, 2) for x in disp_threshold]) + '\n')
        output.write('>Linear Dispersion Method: ' + disp_method + ' (separate linear dispersion threshold for each amino acid group)\n')
    except:
        output.write('>Linear Dispersion Threshold: ' + str(round(disp_threshold, 2)) + '\n')
        output.write('>Linear Dispersion Method: ' + disp_method + ' (a single linear dispersion threshold for all amino acid groups combined)\n')
    output.write('>Composition to Ignore Dispersion: ' + str(round(ign_disp, 2)) + '\n\n')
    
    return output
    
    
def main(args):

    fasta_files, win_size, amino_acids, other_aas, comp_thresholds, disp_threshold, disp_method, ign_disp, uniprots, run_autodetect, impute_value = get_params(args)
    
    random_id = random.randint(0, 10000000)
    
    if uniprots:
        proteome_prots, missing_prots, prots_for_GOanalysis, no_isoforms_prots = crosscheck_user_prots(uniprots, fasta_files)

    if missing_prots:
        print('\nWARNING:\nThe following protein IDs could not be found in your proteome file: ', missing_prots)
        print('Missing proteins could affect LCD enrichment statistics.\n')
        missing_prots_output = open( str(random_id) + '_MissingProteins.txt', 'w')
        missing_prots_output.write('Protein IDs submitted: ' + str(uniprots) + '\n')
        missing_prots_output.write('Proteome FASTA files: ' + str(fasta_files) + '\n')
        missing_prots_output.write('The following protein IDs could not be found in your proteome file:\n')
        missing_prots_output.write('\n'.join(missing_prots))
        missing_prots_output.close()

    if run_autodetect:
        jobfile, lcd_prots_df = LCDcomposer_AutoDetect(random_id, fasta_files, win_size, amino_acids, other_aas, comp_thresholds, disp_threshold, disp_method, ign_disp)

        filtered_proteome_prots = {prot for prot in proteome_prots if prot not in uniprots}     # ENSURES THAT THE "BACKGROUND" SET DOES NOT INCLUDE THE USER UNIPROTS. THIS IS CORRECT CONTINGENCY TABLE DESIGN...EFFECTIVELY TESTS: [ [user_LCD_prots, user_nonLCD_prots], [proteome_LCD_prots_ThatAreNotUserProts, proteome_nonLCD_prots_ThatAreNotUserProts] ]
        line_data = []
        pvals = []
        log2ORs = []
        user_lcdprot_idlists = []
        proteome_lcdprot_idlists = []
        color_palette = []
        imputation_checks = []
        prots_for_GOanalysis = set()    # WILL BE A SET WITH ALL LCD-CONTAINING PROTEINS IN THE USER-SUBMITTED PROTEIN LIST
        for aa in residues:
            lcd_prots = lcd_prots_df[aa]
            user_hits, user_nonhits, proteome_hits, proteome_nonhits, log2OR, pval, user_lcdprots, proteome_lcdprots, imputation_used = test_LCD_enrichment(lcd_prots, uniprots, proteome_prots, missing_prots, filtered_proteome_prots, impute_value)
            pvals.append(pval)
            log2ORs.append(log2OR)
            line_data.append( [user_hits, user_nonhits, proteome_hits, proteome_nonhits] )
            user_lcdprot_idlists.append(user_lcdprots)
            proteome_lcdprot_idlists.append(proteome_lcdprots)
            imputation_checks.append(imputation_used)
            if imputation_used:
                color_palette.append('0.8')
            else:
                color_palette.append('#1f77b4')
            for user_lcdprot in user_lcdprots:
                prots_for_GOanalysis.add(user_lcdprot)

        corrected_pvals = holmsidak_correction(pvals)
        
        prots_for_GOanalysis = [prot for prot in prots_for_GOanalysis if prot in no_isoforms_prots]
        all_submittedprots_with_lcds = [prot for prot in prots_for_GOanalysis]
        user_lcdresults_filename = write_LCDresults_SubmittedProteins(random_id, jobfile, fasta_files, win_size, amino_acids, other_aas, comp_thresholds, disp_threshold, disp_method, ign_disp, all_submittedprots_with_lcds)

        log2ORs = [x if not math.isnan(x) else -5000000 for x in log2ORs]
        log2ORs, line_data, corrected_pvals, uncorrected_pvals, sorted_aas, user_lcdprot_idlists, proteome_lcdprot_idlists, color_palette, imputation_checks = zip(*sorted(zip(log2ORs, line_data, corrected_pvals, pvals, residues, user_lcdprot_idlists, proteome_lcdprot_idlists, color_palette, imputation_checks), reverse=True))
        log2ORs = [x if x != -5000000 else float('nan') for x in log2ORs]

        if not impute_value:
            log2ORs = [x if x != -1000 else 'N/A' for x in log2ORs] #REPLACE DUMMY VALUES WITH "N/A" STRING.
            # USES log2ORs HERE==============
            rows = [[sorted_aas[i]] + line_data[i] + [log2ORs[i], uncorrected_pvals[i], corrected_pvals[i], ', '.join(user_lcdprot_idlists[i]), ', '.join(proteome_lcdprot_idlists[i])] for i in range(len(corrected_pvals))]
        else:
            log2ORs_string = [x if not imputation_checks[i] else 'N/A* [' + str(x) + ' biased est.]' for i, x in enumerate(log2ORs)] #REPLACE DUMMY VALUES WITH "N/A" STRING.
            log2ORs_string_rounded = [x if not imputation_checks[i] else 'N/A* [' + str(round(x, 3)) + ' biased est.]' for i, x in enumerate(log2ORs)] #REPLACE DUMMY VALUES WITH "N/A" STRING.
            # USES log2ORs_string_rounded HERE======
            rows = [[sorted_aas[i]] + line_data[i] + [log2ORs_string_rounded[i], uncorrected_pvals[i], corrected_pvals[i], ', '.join(user_lcdprot_idlists[i]), ', '.join(proteome_lcdprot_idlists[i])] for i in range(len(corrected_pvals))]
        
        header_labels = ['LCD Type', '# of Your Proteins with LCD(s)', '# of Your Proteins without LCD(s)', '# of Background Proteins with LCD(s)', '# of Background Proteins without LCD(s)', 'log2(odds ratio)', 'P-value (Uncorrected)']
        if run_autodetect:
            header_labels += ['Holm-Sidak Corrected P-value']
        stats_filename = write_statsfile(str(random_id), rows, header_labels)
        
        user_lcdprot_idlists = [prots if len(prots) > 0 else {'None'} for prots in user_lcdprot_idlists]
        proteome_lcdprot_idlists = [prots if len(prots) > 0 else {'None'} for prots in proteome_lcdprot_idlists]
        rows = [[sorted_aas[i]] + line_data[i] + [log2ORs[i], uncorrected_pvals[i], corrected_pvals[i]] for i in range(len(corrected_pvals))]
        rows = [rows[i] + [', '.join(user_lcdprot_idlists[i]), ', '.join(proteome_lcdprot_idlists[i])] for i in range(len(rows))]
        
        filtered_aas = [sorted_aas[i] for i in range(len(sorted_aas)) if log2ORs[i] != 'N/A' and not math.isnan(log2ORs[i])]
        filtered_pvals = [float(corrected_pvals[i]) for i in range(len(corrected_pvals)) if log2ORs[i] != 'N/A' and not math.isnan(log2ORs[i])]
        filtered_colorpalette = [color_palette[i] for i in range(len(color_palette)) if log2ORs[i] != 'N/A' and not math.isnan(log2ORs[i])]
        filtered_log2ORs = [log2OR for log2OR in log2ORs if log2OR != 'N/A' and not math.isnan(log2OR)]
        
        # NEEDED TO RE-SORT AFTER FILTERING LISTS...nan VALUES DISRUPT SORTING
        if len(filtered_log2ORs) != 0:
            filtered_log2ORs, filtered_aas, filtered_pvals, filtered_colorpalette = zip(*sorted(zip(filtered_log2ORs, filtered_aas, filtered_pvals, filtered_colorpalette), reverse=True))
            plotting_df = {'log2OR':list(filtered_log2ORs), 'P-value':list(filtered_pvals), 'Amino Acid':list(filtered_aas)}
            plot_filename = barplot(str(random_id), plotting_df, filtered_colorpalette)
        else:
            print('\nWARNING:\nNone of your proteins contain an LCD that passes your search parameters. Without including an imputed value of 1 (optional parameter), an enrichment plot will not be generated for this set of proteins even if the option to generate a plot was selected.\n')
    else:
        if disp_method == 'Original':
            job_file, lcd_prots = LCDcomposer(random_id, fasta_files, win_size, amino_acids, other_aas, comp_thresholds, disp_threshold, disp_method, ign_disp)
        else:
            job_file, lcd_prots = LCDcomposer_NewDispMethod(random_id, fasta_files, win_size, amino_acids, other_aas, comp_thresholds, disp_threshold, disp_method, ign_disp)
        filtered_proteome_prots = {prot for prot in proteome_prots if prot not in uniprots}     # ENSURES THAT THE "BACKGROUND" SET DOES NOT INCLUDE THE USER UNIPROTS. THIS IS CORRECT CONTINGENCY TABLE DESIGN...EFFECTIVELY TESTS: [ [user_LCD_prots, user_nonLCD_prots], [proteome_LCD_prots_ThatAreNotUserProts, proteome_nonLCD_prots_ThatAreNotUserProts] ]

        user_hits, user_nonhits, proteome_hits, proteome_nonhits, log2OR, pval, user_lcdprots, proteome_lcdprots, imputation_used = test_LCD_enrichment(lcd_prots, uniprots, proteome_prots, missing_prots, filtered_proteome_prots, impute_value)

        if log2OR == -1000:
            log2OR = 'N/A'
            
        header_labels = ['LCD Type', '# of Your Proteins with LCD(s)', '# of Your Proteins without LCD(s)', '# of Background Proteins with LCD(s)', '# of Background Proteins without LCD(s)', 'log2(odds ratio)', 'P-value (Uncorrected)']

        user_lcdresults_filename = write_LCDresults_SubmittedProteins(random_id, job_file, fasta_files, win_size, amino_acids, other_aas, comp_thresholds, disp_threshold, disp_method, ign_disp, user_lcdprots)

        if len(user_lcdprots) == 0:
            user_lcdprots = {'None'}
        if len(proteome_lcdprots) == 0:
            proteome_lcdprots = {'None'}
        rows = [['_'.join(amino_acids), user_hits, user_nonhits, proteome_hits, proteome_nonhits, log2OR, pval, ', '.join(user_lcdprots), ', '.join(proteome_lcdprots)]]
        stats_filename = write_statsfile(str(random_id), rows, header_labels)

    if len(fasta_files) > 1:
        print('\nWARNING:\nRunning LCD enrichment statistics while including all protein isoforms is generally not recommended since LCD-containing proteins with many isoforms can skew the statistics. Consider running LCD enrichment tests on a proteome containing only one representative isoform per protein.\n')
    
    print('LCD enrichment analysis complete:', str(datetime.datetime.now()), '\n')


def get_args(arguments):
    parser = argparse.ArgumentParser(description='Test for LCD enrichment among a pre-defined set of proteins.')
    
    parser.add_argument('fasta_file', help="""Your sequence file (in FASTA format).""")
    
    parser.add_argument('uniprots_file', help="""Your file containing UniProt IDs for your proteins of interest (one per line).""")
    
    parser.add_argument('-w', '--window_size', type=str, default='20', 
                        help="""Scanning window size used in the initial search. Default=20aa.
                        However, once a domain is initiated, the domain can be extended indefinitely.
                        Only integers (whole numbers) between 5-10000 are valid window sizes""")
                        
    parser.add_argument('-d', '--dispersion', type=str, default='default',
                        help="""Maximum amino acid dispersion threshold. Default=0.5.
                        
                        This parameter represents the degree of separation/de-mixing of amino acid(s) of interest.
                        All regions with normalized amino acid dispersion below this threshold are discarded. 
                        The value must be a decimal between 0.0 and 1.0 (inclusive).
                        
                        High dispersion values indicate well mixed sequences (e.g. GGGXGGGXGGG, GXGXGXGXGXG).
                        Low dispersion values indicate asymmetric/segregated sequences (e.g. XXGGGGGGGGG, XXXXXGGGGGG).
                        
                        Setting the dispersion threshold value = 0.0 will effectively turn off the dispersion threshold parameter.
                        Setting the dispersion threshold value = 1.0 will result in identification of ONLY perfectly mixed sequences.
                        """)
                        
    parser.add_argument('-c', '--composition', type=str, default='40',
                        help="""Composition threshold for defining X-rich regions, where X represents a particular amino acid or group of amino acids. Default = 40.0 (corresponding to 40% composition).
                        
                        Composition is calculated as the fraction of residues of interest in each window.
                        
                        Value must be between 0-100.
                        """)

    parser.add_argument('-a', '--amino_acids', type=str, default='Auto-detect',
                        help="""Amino acid(s) of interest.
                        
                        Simple amino acid criteria should be an unbroken string consisting of the single-letter abbreviation for a single amino acid or a group of amino acids.
                        e.g.
                        Q       (Q composition and spacing vis-a-vis all other residues)
                        QN      (Composition and spacing of QN residues vis-a-vis all other residues)
                        QNST    (Composition and spacing of QNST residues vis-a-vis all other residues)
                        
                        Complex amino acid criteria allow for specification of different composition thresholds for distinct amino acids or distinct groups of amino acids.
                        Complex criteria use "AND" logic, meaning a sequence must contain each of the amino acid components at or above their respective composition thresholds to be considered.
                        Separate amino acids or groups of amino acids should be separated by an underscore '_'.
                        e.g.
                        Q_P     (must have specified Q composition AND must have specified P composition)
                        QN_ST   (must have specified QN composition AND must have specified ST composition)
                        G_R_Y   (must have specified G composition AND must have specified R composition AND must have specified Y composition)
                        
                        Composition thresholds for each of the separated amino acids or groups are specified by the -c parameter in a similar manner
                        e.g.
                        -a Q_P -c 40_20    (must have 40% Q composition AND must have 20% P composition)
                        -a QN_ST -c 60_25    (must have 60% QN composition AND must have 25% ST composition)
                        """)
                        
    parser.add_argument('-i', '--ignore_dispersion_threshold', type=str, default=None,
                        help="""Composition threshold at which the linear dispersion parameter is ignored when identifying and merging regions that pass the composition threshold.
                        
                        Defaults to half the distance between the user-specified composition threshold and 100%.
                        e.g.
                        A user specified composition threshold = 30% --> ignore_dispersion_threshold = 65%
                        A user specified composition threshold = 40% --> ignore_dispersion_threshold = 70%
                        A user specified composition threshold = 80% --> ignore_dispersion_threshold = 90%
                        
                        However, if the optional -i flag is used, it should be accompanied by a user-specified percentage value.
                        e.g. 
                        -i 60  --> ignore dispersion parameter if composition of the sequence is above 60%
                        -i 75  --> ignore dispersion parameter if composition of the sequence is above 75%
                        -i 100  --> ignore dispersion parameter for all sequences
                        """)

    parser.add_argument('-m', '--dispersion_method', type=str, default='New',
                        help="""Method used to define linear dispersion. Allowable values are "New" or "Original".
                        The "New" method allows you to define a separate linear dispersion threshold for each amino acid or group of amino acids.
                        e.g. for a search with "-a QN_Y" setting distinct groups of amino acids, you could set dispersion thresholds of "-d 0.5_0.7", where the dispersion of Q/N residues must be at least 0.5 but the dispersion of Y residues must be at least 0.7. 
                        
                        The "Original" method is the default behavior of the original LCD-Composer algorithm and only defines a single linear dispersion threshold for all amino acid groups combined.
                        e.g. for the same LCD search ("-a QN_Y"), a single dispersion of "-d 0.6" would apply to all groups combined (Q/N/Y).
                        
                        For LCD searches with only a single amino acid or group of amino acids, the behavior of these two methods is identical.""")
                    
    parser.add_argument('-if', '--isoforms_file', type=str, default=None,
                        help="""FASTA format file containing additional isoform sequences that you wish to analyze.
                        
                        If an isoforms file is specified, LCD-Composer and subsequent statistical analyses of LCD enrichment will include these isoforms in all calculations.
                        Please think carefully about whether it is appropriate to include isoforms in your statistical tests, since LCD-containing proteins with many isoforms can skew the statistics.""")

    parser.add_argument('-r', '--run_autodetect', action='store_true',
                        help="""Run "Auto-detect" mode to identify simple LCDs.
                        
                        This will run a separate LCD search for each of the 20 canonical amino acids using the window size, composition threshold, and linear dispersion threshold specified by the user (or using default values).
                        
                        Note that computation time may be much longer for this mode.""")
                        
    parser.add_argument('-iv', '--impute_value', action='store_true',
                        help="""Impute value of 1 for classes of LCDs for which no LCDs were identified in your proteins of interst.
                        
                        This results in biased estimates of the degree of LCD enrichment in your proteins of interest. However, p-values are still calculated based on the unbiased observation (i.e. using the 0 value, in these instances).""")

    args = parser.parse_args(arguments)
    
    return args

if __name__ == '__main__':
    import sys, argparse, random, datetime, math, statistics, os, mpmath
    import matplotlib.pyplot as plt
    from matplotlib import patches
    import seaborn as sns
    from scipy import stats
    args = get_args(sys.argv[1:])
    main(args)
