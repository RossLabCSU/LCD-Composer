
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


all_aas = 'ACDEFGHIKLMNPQRSTVWY'


def LCDcomposer(random_id, prot_id, seq, win_size, amino_acids, other_aas, comp_thresholds, disp_threshold, ign_disp, output):

    combined_aas = ''.join(amino_acids)

    mins_df = {}
    maxs_df = {}

    #REMOVE STOP CODON FROM C-TERMINUS AND WARN USERS IF A SEQUENCE CONTAINS MULTIPLE STOP CODONS
    if seq.count('*') == 1 and seq[-1] == '*':
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

    #DERIVE DOMAIN BOUNDARIES AND CORRESPONDING DOMAIN SEQUENCES
    domain_boundaries = merge_windows(hit_positions, win_size, seq)
    seqs = [seq[ bounds[0] : bounds[1] ] for bounds in domain_boundaries]

    #TRIM TERMINI UNTIL THEY MATCH THE RESIDUE OF INTEREST (OR ONE OF THE RESIDUES OF INTEREST FOR AA GROUPS)
    trimmed_seqs, trimmed_boundaries = trim_termini(seqs, combined_aas, domain_boundaries)

    #CONVERT TO MORE USER-INTUITIVE PROTEIN BOUNDARY NUMBERING
    trimmed_boundaries = ['('+str(bounds[0]+1) + '-' + str(bounds[1])+')' for bounds in trimmed_boundaries]
    final_comps, final_stddevs = calc_final_comps_stddevs(trimmed_seqs, combined_aas, win_size, mins_df, maxs_df, other_aas)

    contains_hits = False
    if len(trimmed_boundaries) > 0:
        for i in range(len(trimmed_seqs)):
            output.write('\t'.join( [prot_id, trimmed_seqs[i], trimmed_boundaries[i], '_'.join(amino_acids), str(round(final_comps[i], 2)), str(round(final_stddevs[i], 4)), '\t'.join([str(round(trimmed_seqs[i].count(aa) / len(trimmed_seqs[i]) * 100, 2)) for aa in all_aas]) ]) + '\n')
        contains_hits = True
        
    return output, contains_hits
    

def LCDcomposer_NewDispMethod(random_id, prot_id, seq, win_size, amino_acids, other_aas, comp_thresholds, disp_thresholds, ign_disp, output):

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

    #REMOVE STOP CODON FROM C-TERMINUS AND WARN USERS IF A SEQUENCE CONTAINS MULTIPLE STOP CODONS
    if seq.count('*') == 1 and seq[-1] == '*':
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

    #DERIVE DOMAIN BOUNDARIES AND CORRESPONDING DOMAIN SEQUENCES
    domain_boundaries = merge_windows(hit_positions, win_size, seq)
    seqs = [seq[ bounds[0] : bounds[1] ] for bounds in domain_boundaries]

    #TRIM TERMINI UNTIL THEY MATCH THE RESIDUE OF INTEREST (OR ONE OF THE RESIDUES OF INTEREST FOR AA GROUPS)
    trimmed_seqs, trimmed_boundaries = trim_termini(seqs, combined_aas, domain_boundaries)

    #CONVERT TO MORE USER-INTUITIVE PROTEIN BOUNDARY NUMBERING
    trimmed_boundaries = ['('+str(bounds[0]+1) + '-' + str(bounds[1])+')' for bounds in trimmed_boundaries]
    final_comps, final_stddevs = calc_final_comps_stddevs_NewMethod(trimmed_seqs, combined_aas, win_size, mins_df, maxs_df, other_aas, amino_acids, other_aas_groupform)
    
    if len(trimmed_boundaries) > 0:
        for i in range(len(trimmed_seqs)):
            output.write('\t'.join( [prot_id, trimmed_seqs[i], trimmed_boundaries[i], '_'.join(amino_acids), str(round(final_comps[i], 2)), final_stddevs[i], '\t'.join([str(round(trimmed_seqs[i].count(aa) / len(trimmed_seqs[i]) * 100, 2)) for aa in all_aas]) ]) + '\n')
    return output
    

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
    

def merge_windows(hit_positions, win_size, seq):
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

    
#===Composition Plotter Functions===================================================================================

def composition_plotter(prot_id, seq, win_size, amino_acids, comp_threshold, resolution, output_type, figure_length, figure_height, threshold_line, color_palette):

    # PREP DICTIONARY WITH NECESSARY AMINO ACIDS
    df = {}
    if amino_acids == 'Auto-detect':
        for aa in all_aas:
            df[aa] = []
    else:
        for group in amino_acids:
            df[group] = []
    
    seq = check_stopcodons(seq)
    df = calc_composition_plotting(seq, df, amino_acids, win_size)
    plotting_aas = set_plotting_AAs(df, amino_acids, comp_threshold)
    
    # PERFORMS mm-to-inches CONVERSION FOR figure_length AND figure_height WHEN FED IN TO plot_comp_disp
    fig = plot_comp_disp(seq, prot_id, df, plotting_aas, win_size, resolution, output_type, figure_length/25.4, figure_height/25.4, amino_acids, comp_threshold, color_palette, threshold_line)
        
    return fig

        
def calc_composition_plotting(seq, df, amino_acids, win_size):

    for i in range(len(seq) - win_size+1):
        window = seq[i:i+win_size]
        
        #SKIP SEQUENCES SHORTER THAN THE WINDOW SIZE
        if len(window) < win_size:
            continue

        if amino_acids == 'Auto-detect':
            for aa in all_aas:
                df[aa].append( window.count(aa) / len(window) * 100 )
        else:
            for group in amino_acids:
                df[group].append( sum([window.count(aa) / len(window) * 100 for aa in group]) )
                
    return df
    

def check_stopcodons(seq):
    
    #REMOVE STOP CODON FROM C-TERMINUS
    if seq[-1] == '*':
        seq = seq[:-1]
    return seq
    
        
def set_plotting_AAs(df, amino_acids, threshold):
    
    if amino_acids == 'Auto-detect':
        new_aas = ''
        for aa in all_aas:
            if max(df[aa]) >= threshold[0]:
                new_aas += aa
        amino_acids = new_aas
        
    return amino_acids


def plot_comp_disp(seq, prot_id, df, plotting_aas, win_size, resolution, output_type, figure_length, figure_height, amino_acids, comp_threshold, color_palette, threshold_line):

    fig, ax = plt.subplots()
    
    if color_palette == 'Seaborn':
        pal = sns.color_palette('colorblind')
        pal += sns.color_palette('pastel')
    else:
        pal = color_palette
        
    linestyles = ['-']*10 + ['--']*10

    # PLOT DOTTED LINE AT COMPOSITION THRESHOLD
    if threshold_line:
        plt.plot((-1000, len(df[plotting_aas[0]])+1000), (threshold_line, threshold_line), linestyle='--', color='0.8')
        
    index = 0
    for aa in plotting_aas:
        plt.plot([i+(win_size/2)+1 for i in range(len(df[aa]))], [x for x in df[aa]], color=pal[index], linestyle=linestyles[index])
        index += 1
    
    xmargin = len(df[aa]) * 0.015
    plt.xlim(-xmargin, len(df[aa])+win_size+xmargin)
    plt.yticks([x for x in range(0, 120, 20)], labels=[0,20, 40, 60, 80, 100], fontname='Arial', fontsize=16)
    plt.xticks(fontname='Arial', fontsize=16)
    plt.ylabel('AA Composition', fontname='Arial', fontsize=18)
    plt.xlabel('Protein Position', fontname='Arial', fontsize=18)
    
    legend_elements = [Line2D([0], [0], color=pal[i], lw=4, label=plotting_aas[i], linestyle=linestyles[i]) for i in range(len(plotting_aas))]
    
    leg = plt.legend(handles=legend_elements, prop={'size':12, 'family':'Arial'}, loc=2, bbox_to_anchor=(1,1))
    leg.set_title('Amino Acid(s)', prop = {'size':12, 'family':'Arial'})

    fig.set_size_inches(figure_length, figure_height)

    return fig


#===END COMPOSITION PLOTTER FUNCTIONS===========================================================
    

def get_params(args):
    """Gather, define, and validate user parameters."""

    allowed_seq_chars = 'ACDEFGHIKLMNPQRSTVWY*'
    fasta_file = args.fasta_file
    try:
        h = open(fasta_file)
        prot_seq = ''
        line = h.readline().rstrip()
        if not line.startswith('>'):
            print('\nError:\nYour file does not appear to be in FASTA format. Please ensure that you have one header line (starting with ">") followed by one protein sequence.\n')
            exit()
        else:
            for line in h:
                if line.startswith('>'):
                    print('\nError:\nYour file does not appear to be in FASTA format. Please ensure that you have one header line (starting with ">") followed by one protein sequence.\n')
                    exit()
                prot_seq += line.rstrip()
    except:
        print('\nError:\nYour file does not appear to be in FASTA format. Please ensure that you have one header line (starting with ">") followed by one protein sequence.\n')
        exit()
        
    #GET USER-SPECIFIED PARAMETERS=================
    prot_id = args.prot_name
    win_size = args.window_size
    try:
        win_size = int(win_size)
        if win_size < 1:
            exit()
    except:
        print('\nError:\nInvalid window size. The window size must be a positive integer.\n')
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
            print('\nError:\nInvalid amino acid(s) of interest. The "Amino Acid(s) of Interest" can only be composed of the 20 canonical amino acids, with underscores to separate amino acid groups.\n')
            exit()
            
    chars = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ*_'
    other_aas = ''.join([x for x in chars if x not in ''.join(amino_acids)])
    comp_thresholds = args.composition.split('_')
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
    if disp_method == 'New':
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
                    print('\nError:\nwhen using the "New" method for linear dispersion, the number of linear dispersion thresholds must match the number of amino acid groups and composition thresholds you specify, each separated by an underscore (e.g. Amino acids: QN_G_Y, Composition thresholds: 40_15_10, Linear dispersion: 0.5_0.7_0.35).\n')
                    exit()
            except:
                print('\nError:\nThe linear dispersion threshold(s) must be numbers between 0.0-1.0 (inclusive).\n')
                exit()
    elif disp_method == 'Original':
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
        
    if len(amino_acids) != len(comp_thresholds):
        print('\nError:\nInvalid composition thresholds for the specified amino acids or amino acid groups. The number of composition thresholds must match the number of amino acids or amino acid groups specified.\n')
        exit()

    if ign_disp < 0-0.000001 or ign_disp > 100+0.000001:
        print('\nError:\nInvalid ignore dispersion threshold. The ignore dispersion threshold must be a number between 0-100 (inclusive)\n')
        exit()
    
    # GET PLOTTING PARAMETERS=========
    include_plot = False
    if args.plot:
        include_plot = True
        img_resolution = args.img_resolution
        try:
            img_resolution = float(img_resolution)
        except:
            print('\nError:\nInvalid image resolution. The image resolution must be a positive integer. Recommended values are between 600-1200, corresponding to dots per inch (dpi).\n')
            exit()


        img_width = args.img_width
        try:
            img_width = float(img_width)
            if img_width < 50 or img_width > 300:
                exit()
        except:
            print('\nError:\nInvalid image width. The image width must be a number between 50-300 (inclusive).\n')
            exit()
            
        img_height = args.img_height
        try:
            img_height = float(img_height)
            if img_height < 50 or img_height > 300:
                exit()
        except:
            print('\nError:\nInvalid image height. The image height must be a number between 50-300 (inclusive).\n')
            exit()

        img_filetype = args.img_filetype
        allowable_imagetypes = ['eps', 'jpeg', 'jpg', 'pdf', 'pgf', 'png', 'ps', 'raw', 'rgba', 'svg', 'svgz', 'tif', 'tiff']
        if img_filetype not in allowable_imagetypes:
            print('\nError:\nThe image filetype you have specified is not supported. The image filetype must be one of the following:')
            print(allowable_imagetypes, '\n')
            exit()
            
        threshold_line = args.threshold_line
        if threshold_line:
            try:
                threshold_line = float(threshold_line)
            except:
                print('\nError:\nThe y-axis threshold line for plotting must be a number between 1-99 (inclusive).\n')
                exit()
            if threshold_line < 1 or threshold_line > 99:
                print('\nError:\nThe y-axis threshold line for plotting must be a number between 1-99 (inclusive).\n')
                exit()

        plotting_aas = args.plotting_aas
        if plotting_aas == 'default':
            plotting_aas = amino_acids[:]
        else:
            plotting_aas = plotting_aas.split('_')
            plotting_aas = [x.upper() for x in plotting_aas]
            aa_str = ''.join(plotting_aas)
            if 0 in [0 for res in aa_str if res not in allowed_seq_chars[:-1] + '_']:
                print('\nError:\nInvalid plotting amino acids. The plotting amino acids can only be composed of the 20 canonical amino acids, with underscores to separate amino acid groups.\n')
                exit()

        color_palette = args.color_palette
        if color_palette:
            color_palette = color_palette.split('_')
            color_palette = ['#' + color for color in color_palette]
            matches = [re.search(r'^#(?:[0-9a-fA-F]{3}){1,2}$', hex) for hex in color_palette]
            if False in matches:
                error_message += '<br />Error: Invalid color palette. The color palette flag must either not be specified (to use the default color palette from the Seaborn library) or be a series of valid color hex codes with the initial hashtag omitted from each, each separated by an underscore ("_"), and must be of the same length as the number of specified amino acid groups.\n'
            elif len(color_palette) != len(plotting_aas):
                print('\nError:\nInvalid color palette. The color palette must either either not be specified (to use the default color palette from the Seaborn library) or be a series of valid color hex codes with the initial hashtag omitted from each, each separated by an underscore ("_"), and must be of the same length as the number of specified amino acid groups.\n')
                exit()
        else:
            color_palette = sns.color_palette()

    else:
        img_resolution, img_width, img_height, img_filetype, threshold_line, color_palette, plotting_aas = ['']*7

    if run_autodetect:
        amino_acids = 'Auto-detect'
        
    return prot_id, prot_seq, win_size, amino_acids, other_aas, comp_thresholds, disp_threshold, disp_method, ign_disp, run_autodetect, include_plot, img_resolution, img_width, img_height, img_filetype, threshold_line, plotting_aas, color_palette
    

def output_params_header(output, random_id, prot_id, prot_seq, win_size, amino_acids, other_aas, comp_thresholds, disp_threshold, disp_method, ign_disp):
    
    output.write('>*RUNTIME PARAMETERS*\n')
    output.write('>Job ID: ' + str(random_id) + '\n')
    output.write('>Protein Name: ' + prot_id + '\n')
    output.write('>Protein Sequence: ' + prot_seq + '\n')
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
    
    prot_id, prot_seq, win_size, amino_acids, other_aas, comp_thresholds, disp_threshold, disp_method, ign_disp, run_autodetect, include_plot, img_resolution, img_width, img_height, img_filetype, threshold_line, plotting_aas, color_palette = get_params(args)
    random_id = random.randint(0, 10000000)

    job_file = str(random_id) + '_LCDcomposer_RESULTS.tsv'
    output = open(job_file, 'w')
    output_params_header(output, random_id, prot_id, prot_seq, win_size, amino_acids, other_aas, comp_thresholds, disp_threshold, disp_method, ign_disp)
    output.write('\t'.join( ['Protein ID','Domain Sequence','Domain Boundaries','Amino Acid(s) of Interest (i.e. AA or groups of AAs used as the LCD search parameter)','Final Domain Composition','Final Domain Linear Dispersion','\t'.join('ACDEFGHIKLMNPQRSTVWY') ]) + '\n')
    
    # THE "NEW" LINEAR DISPERSION METHOD IS NOT AVAILABLE FOR "Auto-detect" MODE SINCE IT ONLY HAS AN EFFECT IF MULTIPLE AMINO ACID GROUPS ARE SPECIFIED (AUTO-DETECT ALWAYS HAS ONLY ONE AMINO ACID GROUP).
    if run_autodetect:
        lcd_hit_categories = []
        residues = 'ACDEFGHIKLMNPQRSTVWY'
        if disp_method == 'New':
            disp_threshold = disp_threshold[0]
        # RUNS LCD-COMPOSER WITH EACH INDIVIDUAL AA AS THE AA OF INTEREST
        for aa in residues:
            output, contains_hits = LCDcomposer(random_id, prot_id, prot_seq, win_size, aa, residues.replace(aa, ''), comp_thresholds, disp_threshold, ign_disp, output)
            if contains_hits:
                lcd_hit_categories.append(aa)
        plotting_aas = lcd_hit_categories[:]
    else:
        if disp_method == 'Original':
            output, contains_hits = LCDcomposer(random_id, prot_id, prot_seq, win_size, amino_acids, other_aas, comp_thresholds, disp_threshold, ign_disp, output)
        else:
            output = LCDcomposer_NewDispMethod(random_id, prot_id, prot_seq, win_size, amino_acids, other_aas, comp_thresholds, disp_threshold, ign_disp, output)
    
    output.close()

    if include_plot:
        if len(plotting_aas) > 0:
            fig = composition_plotter(random_id, prot_seq, win_size, plotting_aas, comp_thresholds, img_resolution, img_filetype, img_width, img_height, threshold_line, color_palette)
            img_filename = str(random_id) + '_CompositionPlot.' + img_filetype
            plt.savefig(img_filename, bbox_inches='tight', dpi=img_resolution, pil_kwargs={'compression':'tiff_lzw'})
            plt.close()
        else:
            print('\nYour protein does not contain an LCD that passes your search parameters. When this occurs and the "Amino acid(s) of interest" is set to "Auto", a composition plot is not generated.\n')


def get_args(arguments):
    parser = argparse.ArgumentParser(description='Single-protein LCD analysis.')
    
    parser.add_argument('fasta_file', help="""Your sequence file (in FASTA format).""")

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
                    
    parser.add_argument('-r', '--run_autodetect', action='store_true',
                        help="""Run "Auto-detect" mode to identify simple LCDs.
                        
                        This will run a separate LCD search for each of the 20 canonical amino acids using the window size, composition threshold, and linear dispersion threshold specified by the user (or using default values).
                        
                        Note that computation time may be much longer for this mode.""")
                        
    parser.add_argument('-p', '--plot', action='store_true',
                        help="""Include a plot for the given protein. By default, no plot is generated.
                        To include a plot, use the "-p" flag with no additional arguments.""")
                        
    parser.add_argument('-ir', '--img_resolution', type=str, default='600',
                        help="""Image resolution for LCD plot in dots per inch (dpi). Default=600.
                        
                        600dpi should be sufficient for most journals (but also depends on image height and width).
                        
                        NOTE: the image file size will depend on the height, width, and resolution.
                        Extremely high resolution values could result in long computation times and very large files.
                        The maximum recommended resolution (and the highest value usually encountered in journal requirements) is 1200dpi.""")
                        
    parser.add_argument('-iw', '--img_width', type=str, default='150',
                        help="""Image width in mm. Default=150mm (~5.9 inches).
                        
                        Acceptable values are between 50-300.""")
                        
    parser.add_argument('-ih', '--img_height', type=str, default='150',
                        help="""Image height in mm. Default=150mm (~5.9 inches).
                        
                        Acceptable values are between 50-300.""")
                        
    parser.add_argument('-if', '--img_filetype', type=str, default='tif',
                        help="""Image filetype. Default=tif""")
                        
    parser.add_argument('-tl', '--threshold_line', type=str, default=None,
                        help="""Threshold line to draw on the plot. Default=No threshold line.
                        
                        To draw a horizontal threshold line, use this flag followed by a value between 1-99 representing a percent composition threshold.
                        e.g. "-tl 50" will draw a horizontal dotted threshold line that intersects the y-axis at 50%.""")
                        
    parser.add_argument('-cp', '--color_palette', type=str, default=None,
                        help="""Color palette for the LCD plot. Default=Seaborn default color palette.
                        
                        To specify your own color palette, use the -cp flag followed by valid hex codes (without the hashtags) separated by underscores.
                        The number of hex codes in your color palette must match the number of amino acids groups specified.
                        e.g. a color palette for a "-a QN_G_Y" LCD search might be "-cp 1f77b4_ff7f0e_2ca02c".""")
                        
    parser.add_argument('-pn', '--prot_name', type=str, default='No_protein_name_provided',
                        help="""Protein name. This is simply a label to use in the data (no other function).""")

    parser.add_argument('-pa', '--plotting_aas', type=str, default='default',
                        help="""Amino acid groups visualized in LCD plot. By default, this will be the groups of amino acids for which LCDs were discovered in the LCD search.
                        
                        You can also specify custom groups of amino acids using this parameters.
                        The required format is identical that used for composition thresholds (amino acid groups separated by underscores).
                        e.g. -pa QN_G_FWY""")

    args = parser.parse_args(arguments)
    
    return args

if __name__ == '__main__':
    import sys, argparse, random, datetime, math, statistics, os, re, matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    import seaborn as sns
    args = get_args(sys.argv[1:])
    main(args)
