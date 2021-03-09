
import pickle

organisms = ['Scerevisiae', 'Celegans', 'Dmelanogaster', 'Drerio', 'Xlaevis', 'Mmusculus', 'Hsapiens']

def main():

    max_threshold = 25

    goterm_prots_df = {}
    for organism in organisms:
        h = open(organism + '.gaf')
        df = {}
        goterm_prots_df[organism] = {}
        
        for line in h:
            items = line.rstrip().split('\t')
            if len(items) < 10:
                continue
            uniprot = items[1]
            go_id = items[4]
            df[go_id] = df.get(go_id, set())
            df[go_id].add( uniprot )
            goterm_prots_df[organism][go_id] = goterm_prots_df[organism].get(go_id, set())
            goterm_prots_df[organism][go_id].add( uniprot )

        h.close()
    
    pf = open('GOterms_dict_with_ActualAssociatedProtSet_as_Keys.dat', 'wb')
    pickle.dump(goterm_prots_df, pf)
    pf.close()
    
        
if __name__ == '__main__':
    main()