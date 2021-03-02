
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
organisms = ['Scerevisiae', 'Celegans', 'Dmelanogaster', 'Drerio', 'Xlaevis', 'Mmusculus', 'Hsapiens']
proteomes = ['UP000002311_Scerevisiae_NoIsoforms', 'UP000001940_6239_Celegans_NoIsoforms', 'UP000000803_7227_Dmelanogaster_NoIsoforms', 'UP000000437_7955_Drerio_NoIsoforms', 'UP000186698_8355_Xlaevis_NoIsoforms', 'UP000000589_10090_Mmusculus_NoIsoforms', 'UP000005640_9606_Hsapiens_NoIsoforms']
    
def main():

    for i in range(len(organisms)):
        proteome = proteomes[i]
        organism = organisms[i]
        
        master_counts = []
        for aa in amino_acids:
            counts = []
            df = {}
            h = open(proteome + '_' + aa + '_LCD-Composer_RESULTS.tsv')
            for i in range(7):
                h.readline()
            header = h.readline()
            for line in h:
                id, seq, bounds, final_comp, final_disp, *comps = line.rstrip().split('\t')
                junk, uniprot, junk = id.split('|')
                start, end = bounds[1:-1].split('-')
                df[uniprot] = df.get(uniprot, {'LCD Class':[], 'Boundaries':[], 'Sequences':[], 'Lines':[]})
                df[uniprot]['LCD Class'].append(aa)
                df[uniprot]['Boundaries'].append( (int(start), int(end)) )
                df[uniprot]['Sequences'].append(seq)
                df[uniprot]['Lines'].append(','.join( [id, seq, aa, bounds, final_comp, final_disp] + comps) )
            h.close()
            
            for secondary_aa in amino_acids:
                h = open(proteome + '_' + secondary_aa + '_LCD-Composer_RESULTS.tsv')
                
                for i in range(7):
                    h.readline()
                header = h.readline()
                hit_prots = set()
                for line in h:
                    id, seq, bounds, final_comp, final_disp, *comps = line.rstrip().split('\t')
                    junk, uniprot, junk = id.split('|')
                    start, end = bounds[1:-1].split('-')
                    if uniprot in df:
                        for existing_bounds in df[uniprot]['Boundaries']:
                            existing_positions = [x for x in range(existing_bounds[0], existing_bounds[1]+1)]
                            current_lcd_positions = [x for x in range(int(start), int(end)+1)]
                            check_overlap = [False if x not in existing_positions else True for x in current_lcd_positions]    #CHECKS FOR ANY OVERLAP BETWEEN TWO DOMAINS
                            if True not in check_overlap:
                                hit_prots.add(uniprot)
                            else:
                                continue
                h.close()
                                    
                if len(hit_prots) > 0 and aa != secondary_aa:
                    output = open(organism + '_' + aa + '_' + secondary_aa + '_MultiLCD_Proteins.txt', 'w')
                    for prot in hit_prots:
                        output.write(prot + '\n')
                    output.close()
                    
                if secondary_aa == aa:
                    counts.append( 0 )
                else:
                    counts.append( len(hit_prots) )
            master_counts.append( counts )


if __name__ == '__main__':
    main()
