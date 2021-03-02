
import pickle
organisms = ['Scerevisiae', 'Celegans', 'Dmelanogaster', 'Drerio', 'Xlaevis', 'Mmusculus', 'Hsapiens']

def main():

    df = {}
    
    for organism in organisms:
        h = open(organism + '.gaf')
        go_terms = []
        for line in h:
            if line[0] == '!':
                continue
            
            items = line.rstrip().split('\t')
            items = [x for x in items if x != '']
            go_id = [x for x in items if x[0:3] == 'GO:'][0]
            go_terms.append(go_id)

        df[organism] = list(set(go_terms))

    pf = open('All_Organisms_GOids.dat', 'wb')
    pickle.dump(df, pf)
    pf.close()
    

if __name__ == '__main__':
    main()