import sys
dbf = open(sys.argv[1], 'r')
db = []
dbf.readline()
for i in dbf:
        db.append(i.strip())

dfh = open("TCGA-OV_TPM_GeneLevel_geneID.txt", 'r')
rfh = open(sys.argv[1][:-4] + '_TCGA-OV.txt', 'w')

rfh.write(dfh.readline())
for i in dfh:
        line = i.strip().split('\t')
        if line[0] in db:
                db.remove(line[0])
                rfh.write('\t'.join(line) + '\n')
                rfh.flush()
print(db)
print(len(db))
