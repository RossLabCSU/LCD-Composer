
organisms = ['Scerevisiae', 'Celegans', 'Dmelanogaster', 'Drerio', 'Mmusculus', 'Hsapiens']

def main():
    
    exclude_codes = ['ISS', 'ISO', 'ISA', 'ISM', 'IEA']
    for organism in organisms:
        h = open(organism + '.gaf')
        output = open(organism + '_' + '_'.join(exclude_codes) + '_Excluded.gaf', 'w')
        for line in h:
            if line.startswith('!'):
                continue
            items = line.rstrip().split('\t')
            evidence_code = items[6]
            if evidence_code in exclude_codes:
                continue
                
            output.write(line)
        output.close()
            
        h.close()
    

if __name__ == '__main__':
    main()