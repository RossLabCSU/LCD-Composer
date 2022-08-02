
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


def extract_features(query_seq, num_features):
    amino_acids = list('ACDEFGHIKLMNPQRSTVWY')
    comps = [query_seq.count(aa)/len(query_seq) * 100 for aa in amino_acids]
    comps, amino_acids = zip(*sorted(zip(comps, amino_acids), reverse=True))
    thresholds = [x-5 for x in comps[:num_features] if x-5>0.000001]
    residues = list(amino_acids[:len(thresholds)])

    residues, thresholds = check_tied_AAs(amino_acids, comps, num_features, thresholds, residues)

    num_features = len(residues)    # ADJUSTS num_features IF NECESSARY AFTER TIE HANDLING.

    disp_thresholds = []
    for i in range(num_features):
        norm_stddev = calc_dispersion(query_seq, residues[i], 'ACDEFGHIKLMNPQRSTVWY'.replace(residues[i], ''), {}, {})
        disp_threshold = max(norm_stddev - 0.2, 0.0)
        disp_thresholds.append(disp_threshold)
    win_size = max(20, int(len(query_seq)/2))
    
    return '_'.join([str(x) for x in thresholds]), '_'.join(residues), disp_thresholds, win_size
    
    
def check_tied_AAs(amino_acids, comps, num_features, thresholds, residues):

    tied_residues = [amino_acids[i] for i in range(0, len(amino_acids)) if comps[i] == comps[num_features-1]]
    
    # ONLY COMPLETE TIED RESIDUE ANALYSIS IF:
    # 1)THE FINAL RANKED AA FEATURE IS PART OF THE TIE (e.g. residues=['Q', 'N', 'G'] AND tied_residues=['G', 'S'])
    # 2)IF THE TIED RESIDUES ARE NOT EXACTLY EQUAL TO THE TOP num_features NUMBER OF RANKED AAs. (e.g. residues=['Q', 'N', 'G'] AND tied_residues=['Q', 'N', 'G'] AND num_features=3). IN THIS CASE, EVEN THOUGH THE RESIDUES ARE TIED, THE num_features IS LARGE ENOUGH (>=len(residues)) THAT THEY CAN ALL BE CONSIDERED SEPARATE SEARCH FEATURES.
    check = sum([1 for res in tied_residues if res in residues[:num_features]])
    if residues[-1] not in tied_residues or check == len(tied_residues):
        return residues, thresholds

    if len(tied_residues) > 1:
        indices = [amino_acids.index(res) for res in tied_residues]
        tied_thresholds = [comps[index] - 5 for index in indices]
        combined_threshold = sum(tied_thresholds)
        
        # OVERRIDE FINAL FEATURE AND ASSOCIATED COMPOSITION THRESHOLD IF THERE IS A TIE.
        thresholds[-1] = combined_threshold
        residues[-1] = ''.join(tied_residues)
        
        # ENSURES THAT STAND-ALONE RESIDUES THAT HAVE BEEN COMBINED WITH THE TIED RESIDUES DO NOT REMAIN AS INDIVIDUAL RESIDUES.
        thresholds = [threshold for i, threshold in enumerate(thresholds) if residues[i] not in tied_residues]
        residues = [res for res in residues if res not in tied_residues]

    return residues, thresholds
    

def main(args):

    all_aas = 'ACDEFGHIKLMNPQRSTVWY'
    query_seq, win_size, amino_acids, other_aas, comp_thresholds, disp_thresholds, ign_disp, fasta_files = get_params(args)
    combined_aas = ''.join(amino_acids)
    random_id = random.randint(0, 10000000)
    print('Start time:', str(datetime.datetime.now()))
    
    other_aas_groupform = []
    for aa_group in amino_acids:
        other_groupform = all_aas[:]
        for aa in aa_group:
            other_groupform = other_groupform.replace(aa, '')
        other_aas_groupform.append(other_groupform)

    mins_df = {}
    maxs_df = {}
    
    output = open(str(random_id) + '_LCDcomposer_RESULTS.tsv', 'w')
    output_params_header(output, random_id, query_seq, fasta_files, win_size, amino_acids, other_aas, comp_thresholds, disp_thresholds, ign_disp)
    output.write('\t'.join( ['LCD Similarity Rank','Protein Description','UniProt ID (when applicable)','Domain Sequence (Complete LCD)','LCD Similarity Score (Manhattan Distance)','LCD Similarity Score (Euclidean Distance)','Domain Boundaries','Final Domain Composition','Final Domain Linear Dispersion','Best Window Similarity Score, Manhattan Distance','Best Window Sequence, Manhattan Distance','Best Window Similarity Score, Euclidean Distance','Best Window Sequence, Euclidean Distance','\t'.join([aa for aa in all_aas]) ]) + '\n')
    
    proteome_size = get_proteome_size(fasta_files[0])
    if len(fasta_files) > 1:
        proteome_size += get_proteome_size(fasta_files[1])
    
    total_seqs = 0
    perc_range = [int(proteome_size*x/10) for x in range(1, 10)]

    output_lines = []
    hit_lcd_seqs = []
    
    for file_index, file in enumerate(fasta_files):
        h = open(file)

        for id, seq in fasta_parser(h):
        
            if id.count('|') == 2:
                junk, uniprot, junk = id.split('|')
            else:
                uniprot = 'N/A'
                
            #PRINT PROGRESS AT 10% MILESTONES
            if total_seqs in perc_range:
                print(str(round(total_seqs/proteome_size*100)) + '% Complete' + '\t' + str(datetime.datetime.now()))
            
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
            
            if len(trimmed_boundaries) > 0:
                for i in range(len(trimmed_seqs)):
                    output_lines.append( [id, uniprot, trimmed_seqs[i], trimmed_boundaries[i], str(round(final_comps[i], 2)), final_stddevs[i], '\t'.join([str(round(trimmed_seqs[i].count(aa) / len(trimmed_seqs[i]) * 100, 2)) for aa in all_aas]) ])
                    hit_lcd_seqs.append(trimmed_seqs[i])
            else:
                continue
            
        h.close()
    
    output = rank_LCDs(query_seq, hit_lcd_seqs, win_size, output_lines, output)
    
    output.close()
    
    print('End time:', str(datetime.datetime.now()))


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
    
    
def rank_LCDs(query_seq, hit_lcd_seqs, win_size, output_lines, output):
    """Ranks discovered LCDs based on similarity (Manhattan distance) to the query LCD sequence.
    Also scores each LCD by Euclidean distance, but ranking is performed based on Manhattan distance.
    Writes the results to the output file.
    
    Returns:
        output (file handle)
    """
    
    scores = [[], []]
    max_core_scores = [[], []]
    max_core_seqs = [[], []]
    for i, method in enumerate(['Manhattan', 'Euclidean']):
        max_dissimilarity = determine_max_dissimilarity(query_seq, method)
        for lcd_seq in hit_lcd_seqs:
            max_core_score, max_core_seq = determine_max_core(lcd_seq, query_seq, win_size, method)
            max_core_scores[i].append(max_core_score / max_dissimilarity * 100)
            max_core_seqs[i].append(max_core_seq)
            similarity_score = calc_similarity_Distance(query_seq, lcd_seq, method)
            scores[i].append(similarity_score / max_dissimilarity * 100)

    # THIS if CONDITIONAL PREVENTS AN ERROR THAT WOULD OCCUR IF NO LCDs WERE FOUND IN THE LCD SEARCH (RESULTING IN EMPTY LISTS FOR scores AND lines).
    if len(scores[0]) > 0:
        scores_Manhattan, max_core_scores_Manhattan, max_core_seqs_Manhattan, scores_Euclidean, max_core_scores_Euclidean, max_core_seqs_Euclidean, lines = zip(*sorted(zip(scores[0], max_core_scores[0], max_core_seqs[0], scores[1], max_core_scores[1], max_core_seqs[1], output_lines)))
        
        for i in range(len(lines)):
            output.write('\t'.join([str(i+1)] + lines[i][:3] + [str(scores_Manhattan[i]), str(scores_Euclidean[i])] + lines[i][3:6] + [str(max_core_scores_Manhattan[i]), max_core_seqs_Manhattan[i], str(max_core_scores_Euclidean[i]), max_core_seqs_Euclidean[i]] + lines[i][6:]) + '\n')

    return output
    
    
def determine_max_core(lcd_seq, query_seq, win_size, method):

    max_core_score = 0
    max_core_seq = ''
    range_val = max(1, len(lcd_seq)-win_size+1)
    for i in range(range_val):
        window = lcd_seq[i:i+win_size]
        core_score = calc_similarity_Distance(query_seq, window, method)
        if core_score > max_core_score:
            max_core_score = core_score
            max_core_seq = window[:]
            
    return max_core_score, max_core_seq
            
            
def calc_similarity_Distance(string1, string2, method):
    """Function that calculates the total composition similarity score based on one of two metrics.
    Manhattan Distance - This score is calculated by:
    1) calculating the % composition for an amino acid in the query sequence and subject sequence.
    2) calculating the absolute value of the difference in % compositions calculated in step 1.
    3) repeating steps 1-2 for all amino acids found in either the query LCD sequence or the subject LCD sequence.
    4) summing the values from step 3.
    
    Euclidean Distance - This score is determined by:
    1) calculating the % composition for an amino acid in the query sequence and subject sequence.
    2) calculating the square of the difference in % compositions calculated in step 1.
    3) repeating steps 1-2 for all amino acids found in either the query LCD sequence or the subject LCD sequence.
    4) summing the values from step 3.
    5) taking the square-root of the value from step 4.
    
    Returns:
        similarity score (float)
    """
    if method == 'Manhattan':
        similarity_score = sum( [ abs( (string1.count(aa)/len(string1)) - (string2.count(aa)/len(string2)) ) *100 for aa in set(string1+string2) ] )
    else:
        similarity_score = math.sqrt( sum( [ (( (string1.count(aa)/len(string1)) - (string2.count(aa)/len(string2)) ) *100)**2 for aa in 'ACDEFGHIKLMNPQRSTVWY' ] ) )

    return similarity_score

    
def determine_max_dissimilarity(query_seq, method):
    """Maximum possible dissimilarity would occur when a subject sequence is composed of entirely different 
    amino acids than the query sequence. This can be modeled by a situation where the subject sequence 
    is a homopolymer of an amino acid not found in the query sequence.
    e.g. query_sequence = "NNNQYAGGNNNQNQNYANYQYQNNY", subject_sequence = "SSSSSSSSSSSSSSSSSSSSSSSS"

    In rare cases where all 20 amino acids are found in the query sequence, the most dissimilar sequence is 
    a homopolymer of the least common amino acid in the query sequence.
    """
    
    aas = 'ACDEFGHIKLMNPQRSTVWY'
    counts = [query_seq.count(aa) for aa in aas]
    min_index = counts.index( min(counts) )
    leastcommon_aa = aas[min_index]
    dissimilar_seq = leastcommon_aa * len(query_seq)
    max_dissimilarity = calc_similarity_Distance(query_seq, dissimilar_seq, method)
    
    return max_dissimilarity

    
def get_params( args ):
    """Gather, define, and validate user parameters."""
    
    #GET USER-SPECIFIED PARAMETERS=================
    fasta_files = [args.fasta_file]
    isoforms_file = args.isoforms_file
    if isoforms_file:
        fasta_files.append(isoforms_file)
        
    error_message = ''
    query_seq_file = args.query_sequence
    try:
        h = open(query_seq_file)
        query_seqs = []
        for id, seq in fasta_parser(h):
            
            is_invalid = sum([1 for x in set(seq) if x in 'ACDEFGHIKLMNPQRSTVWY']) != len(set(seq))
            if is_invalid:
                print('\nError:\nYour query LCD sequence must be in FASTA format and the sequence can only contain the 20 canonical amino acids (no other symbols are permitted).\n')
                exit()
                
            query_seqs.append(seq)
            
        if len(query_seqs) > 1:
            print('\nError:\nYour query LCD sequence file should only contain one query sequence.\n')
            exit()
        query_seq = query_seqs[0]
    except:
        print('\nError:\nYour query LCD sequence file could not be found. Common causes for this error include (but are not limited to):\n1) The query seqeuence file is not in the same folder/directory as this script,\n2) the query sequence file is misspelled,\n3) there are spaces in the name of the query sequence file, which must be handled correctly in the command (typically surrounded in quotes, but this may be OS-specific)\n')
        exit()

    num_features = args.n_features

    if len(set(query_seq)) < num_features:
        print('\nError:\nThe number of defining compositional features cannot be less than the number of unique amino acids in your query sequence.\n')
        exit()
    else:
        comp_thresholds, amino_acids, disp_thresholds, win_size = extract_features(query_seq, num_features)

    amino_acids = amino_acids.split('_')
    amino_acids = [x.upper() for x in amino_acids]
    chars = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ*_'
    other_aas = ''.join([x for x in chars if x not in ''.join(amino_acids)])
    comp_thresholds = comp_thresholds.split('_')
    comp_thresholds = [float(x)-0.000001 for x in comp_thresholds]    # +/-0.000001 --> ADJUST VALUE FOR DOWNSTREAM FLOATING POINT COMPARISONS
    if not args.ignore_dispersion_threshold:
        ign_disp = sum(comp_thresholds) + ( (100 - sum(comp_thresholds)) / 2 ) - 0.000001
    else:
        ign_disp = args.ignore_dispersion_threshold - 0.000001

    if ign_disp < 0-0.000001 or ign_disp > 100+0.000001:
        print('\nError:\nInvalid ignore dispersion threshold. The ignore dispersion threshold must be a number between 0-100 (inclusive)\n')
        exit()

    return query_seq, win_size, amino_acids, other_aas, comp_thresholds, disp_thresholds, ign_disp, fasta_files
    

def output_params_header(output, random_id, query_seq, fasta_files, win_size, amino_acids, other_aas, comp_thresholds, disp_thresholds, ign_disp):

    output.write('>*RUNTIME PARAMETERS*\n')
    output.write('>Job ID: ' + str(random_id) + '\n')
    output.write('>Query LCD Sequence: ' + query_seq + '\n')
    output.write('>FASTA File(s): ' + ', '.join(fasta_files) + '\n')
    output.write('>Window Size: ' + str(win_size) + '\n')
    output.write('>Amino Acid(s): ' + str(amino_acids) + '\n')
    output.write('>Composition Threshold(s): ' + str([round(x, 1) for x in comp_thresholds]) + '\n')
    output.write('>Linear Dispersion Threshold: ' + str([round(x, 2) for x in disp_thresholds]) + '\n')
    output.write('>Composition to Ignore Dispersion: ' + str(round(ign_disp, 2)) + '\n\n')
    
    return output


def get_args(arguments):
    parser = argparse.ArgumentParser(description='Automated LCD similarity search.\n')
    
    parser.add_argument('fasta_file', help="""Your sequence file (in FASTA format).""")
    
    parser.add_argument('query_sequence', help="""The name of the FASTA file containing your query LCD sequence.""")

    parser.add_argument('-n', '--n_features', type=int, default=2,
                        help="""Number of compositional features to extract from your query LCD sequence and use in the LCD search. Default=2.
                        
                        e.g.
                        An n_features=2 and query LCD sequence of "QQQQPPQPQQPQMPQLQQTQSPPQQQPQQPQPQQ" would result in 
                        an LCD search based on the two most prominent amino acids, Q and P, respectively.
                        
                        The window size is determined based on the length of the query sequence.
                        
                        Separate composition thresholds and linear dispersion thresholds are automatically 
                        calculated for each compositional feature.
                        
                        For detailed descriptions of these calculations, please visit: http://lcd-composer.bmb.colostate.edu:9000/help
                        """)
                        
    parser.add_argument('-i', '--ignore_dispersion_threshold', type=float, default=None,
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
                        
    parser.add_argument('-if', '--isoforms_file', type=str, default=None,
                        help="""FASTA format file containing additional isoform sequences that you wish to analyze.
                        
                        If an isoforms file is specified, LCD-Composer and subsequent statistical analyses of LCD enrichment will include these isoforms in all calculations.
                        Please think carefully about whether it is appropriate to include isoforms in your statistical tests, since LCD-containing proteins with many isoforms can skew the statistics.""")

    args = parser.parse_args(arguments)
    
    return args

if __name__ == '__main__':
    import sys, argparse, random, datetime, math, statistics, os
    args = get_args(sys.argv[1:])
    main(args)