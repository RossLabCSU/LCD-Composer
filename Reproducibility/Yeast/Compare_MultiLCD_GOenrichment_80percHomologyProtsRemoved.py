
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
aa_names = {'A':'Alanine (A)', 'C':'Cysteine (C)', 'D':'Aspartic Acid (D)', 'E':'Glutamic Acid (E)',
        'F':'Phenylalanine (F)', 'G':'Glycine (G)', 'H':'Histidine (H)', 'I':'Isoleucine (I)', 
        'K':'Lysine (K)', 'L':'Leucine (L)', 'M':'Methionine (M)', 'N':'Asparagine (N)',
        'P':'Proline (P)', 'Q':'Glutamine (Q)', 'R':'Arginine (R)', 'S':'Serine (S)',
        'T':'Threonine (T)', 'V':'Valine (V)', 'W':'Tryptophan (W)','Y':'Tyrosine (Y)'}
        
proteome = 'UP000002311_Scerevisiae_NoIsoforms'
organism = 'Scerevisiae'

def main():

    primary_gos = get_primary_GOs()
    secondary_gos = get_secondary_GOs()
    
    primary_gos = pd.DataFrame.from_dict( primary_gos )
    secondary_gos = pd.DataFrame.from_dict( secondary_gos )

    output = open('MultiLCD_Proteins_GOenrichment_80percHomologyProtsRemoved.tsv', 'w')
    output.write('\t'.join(['GO term', 'Primary LCD Class', 'Secondary LCD Class', 'Effect of Subclassing']) + '\n')
    prim_gos = primary_gos[primary_gos['Proteome'] == organism]
    sec_gos = secondary_gos[secondary_gos['Proteome'] == organism]
    
    for aa in amino_acids:
        prim_temp = prim_gos[prim_gos['Primary AA'] == aa]
        sec_temp = sec_gos[sec_gos['Primary AA'] == aa]
        for secondary_aa in amino_acids.replace(aa, ''):
            p_gos = prim_temp[prim_temp['Secondary AA'] == secondary_aa]
            s_gos = sec_temp[sec_temp['Secondary AA'] == secondary_aa]
            
            for go_term in list(p_gos['GO term']):
                if go_term in list(s_gos['GO term']):
                    effect = 'Retained'
                else:
                    effect = 'Lost'
                output.write('\t'.join([go_term, aa, secondary_aa, effect]) + '\n')
                
            for go_term in list(s_gos['GO term']):
                if go_term not in list(p_gos['GO term']):
                    effect = 'New'
                    output.write('\t'.join([go_term, aa, secondary_aa, effect]) + '\n')
                    
    output.close()

    
def get_primary_GOs():
    
    df = {'Proteome':[],
        'GO term':[],
        'Primary AA':[],
        'Secondary AA':[],
        'GO ID':[],
        'Category':[]
        }
        
    for aa in amino_acids:
        for secondary_aa in amino_acids.replace(aa, ''):
            try:
                h = open(aa + '_' + secondary_aa + 'secondary_MultiLCD_Proteins_GO_RESULTS.tsv')
            except:
                continue
            header = h.readline()
            for line in h:
                items = line.rstrip().split('\t')
                if len(items) == 13:
                    go_id, cat, e_or_p, go_desc, ratio_in_list, ratio_in_pop, uncorr_pval, depth, num_prots, bonf_pval, sidak_pval, holm_pval, assoc_prots = items    #lines with 'p' for e_or_p sometimes do not have any proteins in the last cell
                else:
                    go_id, cat, e_or_p, go_desc, ratio_in_list, ratio_in_pop, uncorr_pval, depth, num_prots, bonf_pval, sidak_pval, holm_pval = items
                
                if e_or_p == 'e' and float(sidak_pval) < 0.05:
                    df['Proteome'].append( organism )
                    df['GO term'].append( go_desc )
                    df['GO ID'].append( go_id )
                    df['Primary AA'].append( aa )
                    df['Secondary AA'].append( secondary_aa )
                    df['Category'].append( cat )
                
            h.close()
            
    return df
    
    
def get_secondary_GOs():

    df = {'Proteome':[],
        'GO term':[],
        'Primary AA':[],
        'Secondary AA':[],
        'GO ID':[],
        'Category':[]
        }
        
    for aa in amino_acids:
        for secondary_aa in amino_acids.replace(aa, ''):
            try:
                h = open('Scerevisiae_' + aa + '_' + secondary_aa + '_MultiLCD_GO_RESULTS_80percHomologyProtsRemoved.tsv')
            except:
                continue
            header = h.readline()
            for line in h:
                items = line.rstrip().split('\t')
                if len(items) == 13:
                    go_id, cat, e_or_p, go_desc, ratio_in_list, ratio_in_pop, uncorr_pval, depth, num_prots, bonf_pval, sidak_pval, holm_pval, assoc_prots = items    #lines with 'p' for e_or_p sometimes do not have any proteins in the last cell
                else:
                    go_id, cat, e_or_p, go_desc, ratio_in_list, ratio_in_pop, uncorr_pval, depth, num_prots, bonf_pval, sidak_pval, holm_pval = items
                
                if e_or_p == 'e' and float(sidak_pval) < 0.05:
                    df['Proteome'].append( organism )
                    df['GO term'].append( go_desc )
                    df['GO ID'].append( go_id )
                    df['Primary AA'].append( aa )
                    df['Secondary AA'].append( secondary_aa )
                    df['Category'].append( cat )
                
            h.close()
                
    return df

    
if __name__ == '__main__':
    main()
    
