dfh = open("DNA_methylation_NT84.bg", 'r')
rfh = open("Epimutation_burden_NT84.txt", 'w')
rfh.write('Sample\tEpiBurden\n')
rfh.flush()
header = dfh.readline().strip('\n').split('\t')[87:]

epimut = [0]*84
total_cpg = 0
for i in dfh:
        line = i.strip('\n').split('\t')
        if '.' not in line[3:]: # CpG covered by all samples
                total_cpg += 1

                met = list(map(lambda x: float(x), line[3:]))

                for j in list(range(0,84)):
                        if abs(met[j+84] - met[j]) >= 20:
                                epimut[j] += 1


for i in list(range(0,84)):
        rfh.write(header[i] + '\t' + str(epimut[i] / total_cpg) + '\n')
        rfh.flush()