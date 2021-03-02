
from scipy import stats
import math

amino_acids = 'ACDEFGHIKLMNPQRSTVWY'

def main():
    
    output = open('TableS14_Scerevisiae_MultiLCD_Proteins_GOresults_DegreeOfEnrichment_Analyses.tsv', 'w')
    output.write('\t'.join( ['GO Term', 'Primary LCD Class', 'Secondary LCD Class', 'Primary LCD Ratio in Study', 'Primary LCD Ratio in Population', 'Secondary LCD Ratio in Study', 'Secondary LCD Ratio in Population', 'Primary lnOR', 'Sidak-corrected P-value Primary LCDs', 'Primary Oddsratio', 'Primary OR Lower Confidence Interval', 'Primary OR Upper Confidence Interval', 'Secondary_lnOR', 'Sidak-corrected P-value Primary LCDs', 'Secondary Oddsratio', 'Secondary OR Lower Confidence Interval', 'Secondary OR Upper Confidence Interval', 'Sec-Prim lnOR Difference', 'CIs Overlapping?' ] ) + '\n')
        
    for prim_aa in amino_acids:
        primary_gos = get_prim_GOterm_degrees(prim_aa)
        if len(primary_gos) == 0:
            continue
        for sec_aa in amino_acids.replace(prim_aa, ''):
            
            secondary_gos = get_sec_GOterm_degrees(prim_aa, sec_aa)
            if len(secondary_gos) == 0:
                continue
            
            for go_term in secondary_gos:
                prim_oddsratio, prim_lnOR, prim_upperCI, prim_lowerCI = calc_lnOR(primary_gos[go_term])
                sec_oddsratio, sec_lnOR, sec_upperCI, sec_lowerCI = calc_lnOR(secondary_gos[go_term])
                
                is_CI_overlapping = check_CI_overlap(prim_upperCI, prim_lowerCI, sec_upperCI, sec_lowerCI)
                
                prim_ratio_in_study, prim_ratio_in_pop, prim_pval = primary_gos[go_term]
                sec_ratio_in_study, sec_ratio_in_pop, sec_pval = secondary_gos[go_term]
                
                output.write('\t'.join( [str(x) for x in [ go_term, prim_aa, sec_aa, prim_ratio_in_study, prim_ratio_in_pop, sec_ratio_in_study, sec_ratio_in_pop, prim_lnOR, str(prim_pval), prim_oddsratio, prim_lowerCI, prim_upperCI, sec_lnOR, str(sec_pval), sec_oddsratio, sec_lowerCI, sec_upperCI, sec_lnOR-prim_lnOR, is_CI_overlapping ]] ) + '\n')
            
def check_CI_overlap(prim_upperCI, prim_lowerCI, sec_upperCI, sec_lowerCI):

    overlap = 0
    if prim_upperCI > sec_lowerCI and prim_upperCI < sec_upperCI:
        overlap = 1
        
    if prim_lowerCI > sec_lowerCI and prim_lowerCI < sec_upperCI:
        overlap = 1
    
    return overlap

            
def get_prim_GOterm_degrees(aa):

    try:
        h = open('Scerevisiae_' + aa + '_GO_RESULTS.tsv')
    except:
        return {}
    header = h.readline()
    
    # go_dict is a dictionary with each significantly enriched GO term as the key,
    # and a list containing [ratio_in_study, ratio_in_pop]. NOTE: both are stored as strings,
    # which will need to be parsed later to calculate lnORs.
    go_dict = {}
    
    for line in h:
        items = line.rstrip().split('\t')
        sidak_pval = float(items[10])
        go_term = items[3].replace('"', '')
        ratio_in_study = items[4]
        ratio_in_pop = items[5]
        
        go_dict[go_term] = [ratio_in_study, ratio_in_pop, sidak_pval]
        
    h.close()
        
    return go_dict


def get_sec_GOterm_degrees(aa, sec_aa):

    try:
        h = open(aa + '_' + sec_aa + 'secondary_MultiLCD_Proteins_GO_RESULTS.tsv')
    except:
        return {}
    header = h.readline()
    
    # go_dict is a dictionary with each significantly enriched GO term as the key,
    # and a list containing [ratio_in_study, ratio_in_pop]. NOTE: both are stored as strings,
    # which will need to be parsed later to calculate lnORs.
    go_dict = {}
    
    for line in h:
        items = line.rstrip().split('\t')
        e_or_p = items[2]
        if e_or_p != 'e':
            continue
        sidak_pval = float(items[10])
        if sidak_pval > 0.05+0.00000001:
            continue
        go_term = items[3]
        ratio_in_study = items[4]
        ratio_in_pop = items[5]
        
        go_dict[go_term] = [ratio_in_study, ratio_in_pop, sidak_pval]
        
    h.close()
        
    return go_dict
    
def calc_lnOR(ratio_list):

    study_hits, study_total = ratio_list[0].split('/')
    pop_hits, pop_total = ratio_list[1].split('/')

    study_hits, study_total, pop_hits, pop_total = int(study_hits), int(study_total), int(pop_hits), int(pop_total)
    contingency = [ [study_hits, study_total], [pop_hits, pop_total] ]
    
    oddsratio, pval = stats.fisher_exact(contingency, alternative='two-sided')
    lnOR = math.log(oddsratio)
    
    upper_CI = math.exp( lnOR + 1.96 * math.sqrt(1/study_hits + 1/study_total + 1/pop_hits + 1/pop_total) )
    lower_CI = math.exp( lnOR - 1.96 * math.sqrt(1/study_hits + 1/study_total + 1/pop_hits + 1/pop_total) )

    return oddsratio, lnOR, upper_CI, lower_CI
    


if __name__ == '__main__':
    main()
