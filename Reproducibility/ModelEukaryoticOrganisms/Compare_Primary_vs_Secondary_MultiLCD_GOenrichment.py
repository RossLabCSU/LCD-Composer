
import pandas as pd
import pickle
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
aa_names = {'A':'Alanine (A)', 'C':'Cysteine (C)', 'D':'Aspartic Acid (D)', 'E':'Glutamic Acid (E)',
        'F':'Phenylalanine (F)', 'G':'Glycine (G)', 'H':'Histidine (H)', 'I':'Isoleucine (I)', 
        'K':'Lysine (K)', 'L':'Leucine (L)', 'M':'Methionine (M)', 'N':'Asparagine (N)',
        'P':'Proline (P)', 'Q':'Glutamine (Q)', 'R':'Arginine (R)', 'S':'Serine (S)',
        'T':'Threonine (T)', 'V':'Valine (V)', 'W':'Tryptophan (W)','Y':'Tyrosine (Y)'}

organisms = ['Scerevisiae', 'Celegans', 'Dmelanogaster', 'Drerio', 'Xlaevis', 'Mmusculus', 'Hsapiens']
        
def main():

    primary_gos = get_primary_GOs()
    secondary_gos = get_secondary_GOs()
    
    primary_gos = pd.DataFrame.from_dict( primary_gos )
    secondary_gos = pd.DataFrame.from_dict( secondary_gos )

    output = open('Table S11 - ModelEukaryoticOrganisms_MultiLCD_subclassification_GOtermEffect.tsv', 'w')
    output.write('\t'.join(['GO term', 'Primary LCD Class', 'Secondary LCD Class', 'Effect of Subclassing']) + '\n')
    
    for organism in organisms:
        prim_gos = primary_gos[primary_gos['Proteome'] == organism]
        sec_gos = secondary_gos[secondary_gos['Proteome'] == organism]
        
        for aa in amino_acids:
            p_gos = prim_gos[prim_gos['Amino Acid'] == aa]
            s_gos = sec_gos[sec_gos['Primary AA'] == aa]
            
            for go_term in list(p_gos['GO term']):
                maintained_secondary_aas = []
                effect = 'Lost' #Initializes effect variable...if the go term is not in any of the secondary lists, then this will output as "Lost"
                for sec_aa in amino_acids.replace(aa, ''):
                    final_s_gos = s_gos[s_gos['Secondary AA'] == sec_aa]
                    
                    if go_term in list(final_s_gos['GO term']):
                        effect = 'Retained'
                        maintained_secondary_aas.append( sec_aa )

                if effect == 'Retained':
                    output.write('\t'.join([organism[0]+'. '+organism[1:], go_term, aa, '_'.join(maintained_secondary_aas), effect]) + '\n')
                else:
                    output.write('\t'.join([organism[0]+'. '+organism[1:], go_term, aa, 'None', effect]) + '\n')
                        
            for sec_aa in amino_acids.replace(aa, ''):
                final_s_gos = s_gos[s_gos['Secondary AA'] == sec_aa]
                
                for go_term in list(final_s_gos['GO term']):
                    if go_term not in list(p_gos['GO term']):
                        effect = 'New'
                        output.write('\t'.join([organism[0]+'. '+organism[1:], go_term, aa, sec_aa, effect]) + '\n')
                    
    output.close()

    
def get_primary_GOs():
    
    df = {'Proteome':[],
        'GO term':[],
        'Amino Acid':[],
        'GO ID':[],
        'Category':[]
        }
        
    for aa in amino_acids:
        for organism in organisms:
            all_go_hits = []
            try:
                h = open(organism + '_' + aa + '_GO_RESULTS.tsv')
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
                    df['Amino Acid'].append( aa )
                    df['Category'].append( cat )
            
            h.close()
            
    pf = open('Primary_MultiLCD_GOenrichment_data.dat', 'wb')
    pickle.dump(df, pf)
    pf.close()
            
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
            for organism in organisms:
                try:
                    # h = open(organism + '_' + aa + '-primary_' + secondary_aa + '-secondary_GO_RESULTS.tsv')
                    h = open(organism + '_' + aa + '_' + secondary_aa + '_MultiLCD_Proteins_GO_RESULTS.tsv')
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
                
    pf = open('Secondary_MultiLCD_Subclass_GOenrichment_data.dat', 'wb')
    pickle.dump(df, pf)
    pf.close()
                
    return df

    
if __name__ == '__main__':
    main()
    
