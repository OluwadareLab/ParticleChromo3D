import Helper

import numpy as np
from scipy.spatial import distance
from scipy import stats

#matrixFilePtr = './data/input-and-models/Input/Synthetic/'
matrixFilePtr = './utility/GSM1173492_Th1_cell1/1mb/iced/Iced_chrX__1mb.txt'

#pdbFilePtr = './output/chr1_500kb_matrix_1589758743647.pdb'
#pdbFilePtr = '../OthersResults/code/3DMax-1.0/examples/input/Matrix/input-and-models/Models/Synthetic/Chromosome3D/'
#pdbFilePtr = '../../3dMax/3DMax-1.0/examples/output/Matrix/Synthetic/'
pdbFilePtr = './output/scale/Th1_Cell1_chrX_model.pdb'

resultsPtr = '../OthersResults/Results/SCL_Results_TH1_cell1_chrmX.txt'

tmp = [ 'chainDres5_Matrix_noise000',
		'chainDres5_Matrix_noise030',
		'chainDres5_Matrix_noise050',
		'chainDres5_Matrix_noise070',
		'chainDres5_Matrix_noise100',
		'chainDres5_Matrix_noise130',
		'chainDres5_Matrix_noise150',
		'chainDres5_Matrix_noise170',
		'chainDres5_Matrix_noise200',
		'chainDres5_Matrix_noise250',
		'chainDres5_Matrix_noise300',
		'chainDres5_Matrix_noise400',]

results = []
for i in range(len(tmp)):
	#tmpMatrixPtr = matrixFilePtr + 'chr' + str(i) + '_500kb_matrix.txt'
	#tempPdbPtr = pdbFilePtr + 'chr' + str(i) + '_500kb_rank_reduced.pdb'

	#tmpMatrixPtr = matrixFilePtr + tmp[i] + '.txt'
	#tempPdbPtr = pdbFilePtr + tmp[i] + '_model.pdb'

	tmpMatrixPtr = matrixFilePtr
	tempPdbPtr = pdbFilePtr

	contact, pointMap, zeroInd = Helper.Read_Data(tmpMatrixPtr,1000,3.1)
	xyz = Helper.Read_PDB(tempPdbPtr)
	#xyz, trash = Helper.Read_Data_List(tempPdbPtr)

	print(xyz.shape)

	calcDist = distance.pdist(xyz)
	print(calcDist.shape)
	if zeroInd is not None:
		calcDist = np.delete(calcDist, zeroInd)
	print(calcDist.shape)
	print(contact[:,3].shape)

	pers = stats.pearsonr(calcDist, contact[:,3])
	spear = stats.spearmanr(calcDist, contact[:,3])
	error = np.sqrt( (1/len(pointMap)) * np.sum( (calcDist-contact[:,3])**2 ) )

	print(	' Pearson: ' + str(pers[0]) + 
			' Spearmen: ' + str(spear[0]) +
			' error: ' + str(error))

	asd=asdsad

	results.append((pers[0], spear[0], error))

Helper.Write_List(results, resultsPtr)