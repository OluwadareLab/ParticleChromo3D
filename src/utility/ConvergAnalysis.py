import Helper

import numpy as np
import pandas as pd

from pathlib import Path

def Write_Avg_IF(inFilePtr, outFilePtr):
	avgIF = []
	for k in range(1, 2):
		print('chr ' + (str)(k))
		#tempInPtr = inFilePtr + 'chr' + (str)(k) + '_1mb_matrix.txt'
		tempInPtr = './data/input-and-models/Input/Synthetic/chainDres5_Matrix_noise000.txt'
		#tempInPtr = inFilePtr + 'chr' + (str)(k) + '_500kb_matrix.txt'
		tmpMatrix, tmp = Helper.Read_Matrix_To_List(tempInPtr)
		avgIF.append(np.mean(tmpMatrix[:,2]))

	Helper.Write_List(avgIF, outFilePtr+'avgIFSynth.txt')

def Write_High_Converg(inFilePtr, outFilePtr):
	lineSpace = [-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,2,2.1,2.2,2.3,2.4,2.5,2.6]

	highConvergeInd = []
	highConvergeVal = []
	for k in range(2, 23):
		print('chr ' + (str)(k))
		tempInPtr = inFilePtr + 'chr' + (str)(k) + '_1mb_matrix' + (str)(k) + '_converg.txt'
		tmpList = Helper.Read_Data_List(tempInPtr)
		highIndex = np.argmax(tmpList)

		highConvergeVal.append(tmpList[highIndex])
		highConvergeInd.append(lineSpace[highIndex])

	Helper.Write_List(highConvergeVal, outFilePtr+'highVal.txt')
	Helper.Write_List(highConvergeInd, outFilePtr+'highInd.txt')

def Write_Bead_Count(inFilePtr, outFilePtr):
	count = []
	for k in range(1, 23):
		print('chr ' + (str)(k))
		#tempInPtr = inFilePtr + 'chr' + (str)(k) + '_1mb_matrix.txt'
		#tempInPtr = './data/input-and-models/Input/Synthetic/chainDres5_Matrix_noise000.txt'
		tempInPtr = inFilePtr + 'chr' + (str)(k) + '_500kb_matrix.txt'
		tmpMatrix, points, tmp = Helper.Read_Data(tempInPtr, 0.4)
		count.append(len(points))

	Helper.Write_List(count, outFilePtr+'beadCount500.txt')

def Write_File_Size(inFilePtr, outFilePtr):
	sizeList = []
	for k in range(1, 23):
		print('chr ' + (str)(k))
		tempInPtr = inFilePtr + 'chr' + (str)(k) + '_1mb_matrix.txt'
		#tempInPtr = './data/input-and-models/Input/Synthetic/chainDres5_Matrix_noise000.txt'
		size = Path(tempInPtr).stat().st_size
		sizeList.append(size)

	Helper.Write_List(sizeList, outFilePtr+'fileSize.txt')

def Write_High_Range(inFilePtr, outFilePtr):
	spearConvVal = []
	costConvVal = []
	for k in range(1, 23):
		print('chr ' + (str)(k))
		tempInPtr = inFilePtr + 'chr' + (str)(k) + '_1mb_matrix' + (str)(k) + '_converg.txt'
		tmpList = np.genfromtxt(tempInPtr, delimiter=' ')

		spearIndex = np.argmax(tmpList[:,1])
		costIndex = np.argmin(tmpList[:,2])

		spearConvVal.append(tmpList[spearIndex])
		costConvVal.append(tmpList[costIndex])

		print(spearConvVal)
		print(costConvVal)

	Helper.Write_List(spearConvVal, outFilePtr+'spearConvVal.txt')
	Helper.Write_List(costConvVal, outFilePtr+'costConvVal.txt')

def Write_Avg_IF_Single(inFilePtr, outFilePtr):
	avgIF = []
	for k in range(1, 2):
		print('chr ' + (str)(k))
		#tempInPtr = './output/Full_HiC/Analysis/list/iced_chr1_HiC_1mb.txt'
		tempInPtr = './output/Full_HiC/Analysis/list/iced_nomulti_chr1_HiC_1mb.txt'
		df = pd.read_table(tempInPtr, delimiter=' ', dtype={"bead1": float, "bead2": float, "if": float})
		beadList = df['bead1'].unique()

		for bead in beadList:
			# Get specific chromsome contacts only
			contacts1 = df.loc[df['bead1'] == bead]
			contacts2 = df.loc[df['bead2'] == bead]
			contacts = pd.concat([contacts1, contacts2])

			mean = contacts[['if']].values.mean()

			avgIF.append((bead, mean))

	Helper.Write_List(avgIF, outFilePtr+'iced_nomulti_chr1_avgIFBead.txt')

inMatrixFilePtr = './data/input-and-models/Input/HiC/'
#inDataFilePtr = './output/Full_HiC/'
inDataFilePtr = './output/Matrix/conv/scale/'

outFilePtr = './output/Full_HiC/Analysis/'

#Write_Avg_IF(inMatrixFilePtr, outFilePtr)
#Write_Bead_Count(inMatrixFilePtr, outFilePtr)
#Write_File_Size(inMatrixFilePtr, outFilePtr)
#Write_High_Converg(inDataFilePtr, outFilePtr)
Write_Avg_IF_Single(inDataFilePtr, outFilePtr)