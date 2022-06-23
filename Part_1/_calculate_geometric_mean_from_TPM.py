import sys
dfh = open(sys.argv[1], 'r')
rfh = open('GM_' + sys.argv[1], 'w') # Geometric Mean
line = dfh.readline().strip('\n').split('\t')
rfh.write('\t'.join(line) + '\n')
gm = list(0 for a in range(len(line[1:])))
c = 0
for i in dfh:
        c += 1
        line = i.strip('\n').split('\t')
        rfh.write('\t'.join(line) + '\n')
        if c == 1:
                gm = list(map(lambda x: float(x)+1, line[1:]))
        else:
                gm_next = list(map(lambda x: (float(line[1:][x])+1)*(gm[x]), range(len(line[1:]))))
                gm = gm_next
rfh.write('EMT_index\t' + '\t'.join(list(map(lambda x: str(x**(1.0/c)), gm))) + '\n')
