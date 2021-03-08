
amino_acids = 'ACDEFGHIKLMNPQRSTVY' # W IS EXCLUDED BECAUSE THERE WERE NO PRIMARY LCDs WITH W COMPOSITION >= 40%

def main():
    threshold = 80
    output = open('HighHomologyProts_' + str(threshold) + 'percentThreshold_MultiLCD_Prots.tsv', 'w')
    output.write('Primary LCD Class\tSecondary AA\tProtein\n')
    threshold = threshold - 0.0000001
    for aa in amino_acids:
        for res2 in amino_acids.replace(aa, ''):
            try:
                h = open(aa + '_' + res2 + '_MultiLCD_Prots.pim.pim')
            except:
                continue
            for i in range(6):
                h.readline()
            
            prots = []
            linecount = 2
            for line in h:
                items = line.rstrip().split(' ')
                items = [x for x in items if x != '']
                prot = items[1]
                prots.append(prot)
                toggle = 0
                for seq_ident in items[2:linecount] + items[linecount+1:]:
                    if float(seq_ident) >= threshold and toggle == 0:
                        output.write(aa + '\t' + res2 + '\t' + prot + '\n')
                        toggle = 1
                linecount += 1
                        
            h.close()
    output.close()

if __name__ == '__main__':
    main()