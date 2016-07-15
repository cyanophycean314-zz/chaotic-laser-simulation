#Trim files

rates = [100,250,500,1000,2000,3200,5000,10000,20000,30000,40000,50000,100000]
filelist = [str(x) + "com" for x in rates] + [str(x) + "comNT" for x in rates]
linestokeep = 10000000

for filename in filelist:
	fin = open(filename + 'vt.out','r')
	fout = open(filename + 'tr_vt.out','w')
	fout.write(fin.readline())
	lines = 0
	while lines < linestokeep:
		fout.write(fin.readline())
		lines += 1