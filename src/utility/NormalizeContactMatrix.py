from iced import filter, normalization

import numpy as np

filePtr = './GSM1173492_Th1_ensemble/50kb/'

for i in range(1,11):
	counts = np.genfromtxt(filePtr+'chr'+str(i)+'_50kb.txt', delimiter=' ')

	counts = filter.filter_low_counts(counts, percentage=0.04)
	normed = normalization.ICE_normalization(counts)

	np.savetxt(filePtr+'Iced_chr'+str(i)+'_'+'_50kb.txt', normed, delimiter=',')
