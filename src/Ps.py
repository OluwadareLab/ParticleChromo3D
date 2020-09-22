from Swarm import Swarm
import Helper

import sys
import copy
import time

import numpy as np
from scipy import stats

from multiprocessing.dummy import Pool
from multiprocessing import cpu_count
from itertools import repeat

def Print_Stats(swarm, contact, pointCount, j, i, outFilePtr):
	pers = stats.pearsonr(swarm.gBest[2], contact[:,3])
	spear = stats.spearmanr(swarm.gBest[2], contact[:,3])
	spearIF = stats.spearmanr(swarm.gBest[2], contact[:,2])

	error = np.sqrt( (1/pointCount) * np.sum( (swarm.gBest[2]-contact[:,3])**2 ) )

	print('id: ' + str(swarm.id) + 
		' itt: ' + str(i) + 
		' Cost: ' + str(swarm.gBest[1]) + 
		' Pearson: ' + str(pers[0]) + 
		' Spearmen: ' + str(spear[0]) +
		' IFSpear: ' + str(spearIF[0]) +
		' error: ' + str(error))

	Write_Stats(swarm, contact, j, outFilePtr)

def Write_Stats(swarm, contact, j, outFilePtr):
	#Helper.Write_List(swarm.gBest[2], outFilePtr + 'testAct' + str(j) + '.txt')
	#Helper.Write_List(contact[:,3], outFilePtr + 'testWish' + str(j) + '.txt')
	#Helper.Write_List(swarm.gBest[0], outFilePtr + 'testPosEnd' + str(j) + '.txt')
	Helper.Write_Output(outFilePtr+str(j), swarm.gBest[0])

def One_Move(ittCount, swarm, contact, pointCount, threshold, j, outFilePtr):
	saveGBestCost = float('inf')
	totTime = 0

	start = time.time()
	for i in range(ittCount):
		if (i%1000 == 0) and (swarm.gBest is not None):
			timeSinceUpdate = time.time()-start
			totTime += timeSinceUpdate
			print(timeSinceUpdate)
			start = time.time()
			print(i, j)
			error = np.sqrt( (1/pointCount) * np.sum( (swarm.gBest[2]-contact[:,3])**2 ) )
			Print_Stats(swarm, contact, pointCount, j, i, outFilePtr)
			if (np.abs(saveGBestCost-error)) >= threshold:
				saveGBestCost = error
			else:
				return i, totTime

		operation(i, swarm)
		
	timeSinceUpdate = time.time()-start
	totTime += timeSinceUpdate
	return i, totTime

def operation(i, swarm):
	swarm.Calc_Vel(ittCount,i)
	swarm.Update_Pos(i)
	swarm.Cost()

def Optimize(maxRange, inFilePtr, outFilePtr):
	contact, points, zeroInd = Helper.Read_Data(inFilePtr,maxRange, 3.1)
	swarm = Swarm(contact, len(points), randVal=randRange, swarmCount=swarmCount, zeroInd=zeroInd)
	Helper.Write_List(swarm.pos[0], './output/testPosStart' + str(maxRange) + '.txt')

	ittFin, totTime = One_Move(ittCount, swarm, contact, len(points), threshold, maxRange, outFilePtr)
	return (stats.pearsonr(swarm.gBest[2], contact[:,3])[0], 
					stats.spearmanr(swarm.gBest[2], contact[:,3])[0], 
					np.sqrt( (1/len(points)) * np.sum( (swarm.gBest[2]-contact[:,3])**2 ) ),
					ittFin,
					totTime, 
					maxRange, swarm.id)

def Par_Choice(lineSpace, inFilePtr, outFilePtr, k):

	bestSwarm = None
	if lineSpace is None:
		bestSwarm = Optimize(None, inFilePtr, outFilePtr)
	elif len(lineSpace) > 1:
		convStore = []
		pool = Pool(processes=(PROC_COUNT))
		swarms = pool.starmap(Optimize, zip(range(lineSpace[0],lineSpace[1],5000), repeat(inFilePtr), repeat(outFilePtr)))

		pool.close()
		pool.join()

		swarms = sorted(swarms, key=lambda x: x[1])

		for swarm in swarms:
			print(str(swarm[-1]) + ' ' + str(swarm[1]))
			convStore.append(swarm)
			if (bestSwarm is None) or (swarm[1] > bestSwarm[1]):
				bestSwarm = swarm
	else:
		bestSwarm = Optimize(lineSpace[0], inFilePtr, outFilePtr)

	print(bestSwarm)

	#Helper.Write_List(convStore, outFilePtr + (str)(k) + '_converg.txt')
	return bestSwarm

def Full_List(countMin, countMax):
	#lineSpace = [-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1,1.1,1.2,1.3,1.4,1.5,1.6,2,2.1,2.2,2.3,2.4,2.5,2.6,3,3.1,3.2,3.3,3.4,3.5,3.6]
	#lineSpace = [0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4]
	#lineSpace = [0.2,0.3,0.4,0.5,0.6,0.7]
	#lineSpace = [3.1]
	#lineSpace = None
	lineSpace = [5000, 40000]
	#lineSpace = [5,5.1,5.2,5.3,5.4,5.5,5.6,5.7,5.8,5.9,6,6.1,6.2,6.3,6.4,6.5,6.6,4,4.1,4.2,4.3,4.4,4.5,4.6]

	tmp = [ 'chainDres5_Matrix_noise000.txt',
			'chainDres5_Matrix_noise030.txt',
			'chainDres5_Matrix_noise050.txt',
			'chainDres5_Matrix_noise070.txt',
			'chainDres5_Matrix_noise100.txt',
			'chainDres5_Matrix_noise130.txt',
			'chainDres5_Matrix_noise150.txt',
			'chainDres5_Matrix_noise170.txt',
			'chainDres5_Matrix_noise200.txt',
			'chainDres5_Matrix_noise250.txt',
			'chainDres5_Matrix_noise300.txt',
			'chainDres5_Matrix_noise400.txt',]

	inFilePtr = './data/input-and-models/Input/HiC/'
	#inFilePtr = './data/input-and-models/Input/Synthetic/'
	#inFilePtr = './utility/GSM1173492_Th1_ensemble/50kb/iced/'

	outFilePtr = './output/Matrix/'
	#outFilePtr = './output/GSM1173492_Th1_ensemble/'
	

	convStore = []
	for k in range(countMin, countMax):
		#for j in range(2,len(tmp)):
		print('chr ' + (str)(k))
		tempInPtr = inFilePtr + 'chr' + (str)(k) + '_500kb_matrix.txt'
		tempOutPtr = outFilePtr + 'chr' + (str)(k) + '_500kb_matrix'
		#tempInPtr = inFilePtr + tmp[j]
		#tempOutPtr = outFilePtr + tmp[j]
		#tempInPtr = inFilePtr + 'Iced_chr' + (str)(k) + '__50kb.txt'
		#tempOutPtr = outFilePtr + 'chr' + (str)(k) + '_50kb_matrix'
		#tempInPtr = inFilePtr + 'Iced_chrX__1mb.txt'
		#tempOutPtr = outFilePtr + 'chrX_1mb_matrix'
		print(tempInPtr)
		convStore.append(Par_Choice(lineSpace, tempInPtr, tempOutPtr, k))
		print(convStore)

	Helper.Write_List(convStore, outFilePtr + 'par_50kb_converg.txt')


sys.setrecursionlimit(10000)
PROC_COUNT = cpu_count()

'''#mkl.set_num_threads(56)
mkl.MKL_DYNAMIC=False
mkl.set_num_threads(16)
print(mkl.get_max_threads())'''

swarmCount = 20
ittCount = 20000
threshold = 0.1
randRange = 1

Full_List(1,2)

'''
# 	 Convergence number Selection Times
#	 Multithreading
#	 Grouping0.3
#		

	5 heuristic methods
	simulated
	lcl method
	chromosome3d method
	name
	chimera visualaization
	
	close gap
	Same scale pdb	- 	DONE
	Standardized distance
	finding optimal convergence factor and random values
	implement matrix/tupple switch

	time compare
	RMSE scale
	PSO distance value range

''' 
