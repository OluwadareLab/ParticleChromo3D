from Swarm import Swarm
import Helper

import sys
import copy
import time
import argparse

import numpy as np
from scipy import stats

from multiprocessing.dummy import Pool
from multiprocessing import cpu_count
from itertools import repeat

# Prints statistics of the current swarm
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
	Helper.Write_Output(outFilePtr+str(j), swarm.gBest[0])

# Performs one operation and prints statistics of current swarm
def One_Move(ittCount, swarm, contact, pointCount, threshold, j, outFilePtr):
	saveGBestCost = float('inf')
	totTime = 0

	start = time.time()
	for i in range(ittCount):
		if (i%1000 == 0) and (swarm.gBest is not None):
			timeSinceUpdate = time.time()-start
			totTime += timeSinceUpdate

			start = time.time()
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

# Performs a single PSO pass: Velocity calculation, update position, get new cost
def operation(i, swarm):
	swarm.Calc_Vel(ittCount,i)
	swarm.Update_Pos(i)
	swarm.Cost()

# Optimizes single swarm
def Optimize(maxRange, inFilePtr, outFilePtr, convFact = 3.1):
	contact, points, zeroInd = Helper.Read_Data(inFilePtr,maxRange, convFact)
	swarm = Swarm(contact, len(points), randVal=randRange, swarmCount=swarmCount, zeroInd=zeroInd)

	ittFin, totTime = One_Move(ittCount, swarm, contact, len(points), threshold, maxRange, outFilePtr)
	return (stats.pearsonr(swarm.gBest[2], contact[:,3])[0], 
					stats.spearmanr(swarm.gBest[2], contact[:,3])[0], 
					np.sqrt( (1/len(points)) * np.sum( (swarm.gBest[2]-contact[:,3])**2 ) ),
					ittFin,
					totTime, 
					maxRange, swarm.id)

# Runs in paralel if passed multiple rangeSpace
def Par_Choice(rangeSpace, inFilePtr, outFilePtr):

	bestSwarm = None
	if rangeSpace is None:
		bestSwarm = Optimize(None, inFilePtr, outFilePtr)
	elif len(rangeSpace) > 1:
		convStore = []
		pool = Pool(processes=(PROC_COUNT))
		swarms = pool.starmap(Optimize, zip(range(rangeSpace[0],rangeSpace[1],5000), repeat(inFilePtr), repeat(outFilePtr)))

		pool.close()
		pool.join()

		swarms = sorted(swarms, key=lambda x: x[1])

		for swarm in swarms:
			print(str(swarm[-1]) + ' ' + str(swarm[1]))
			convStore.append(swarm)
			if (bestSwarm is None) or (swarm[1] > bestSwarm[1]):
				bestSwarm = swarm
	else:
		bestSwarm = Optimize(rangeSpace[0], inFilePtr, outFilePtr)

	print(bestSwarm)

	return bestSwarm

def Full_List(rangeSpace, inputFilePtr, outFilePtr):
	convStore = []

	print(tempInPtr)
	convStore.append(Par_Choice(rangeSpace, inputFilePtr, outFilePtr))
	print(convStore)

	Helper.Write_List(convStore, outFilePtr)


sys.setrecursionlimit(10000)
PROC_COUNT = cpu_count()

inFilePtr = '../input-and-models/Input/HiC/'
outFilePtr = './chr.pdb"'

tempInPtr = inFilePtr + 'chr' + str(12) + '_1mb_matrix.txt'
tempOutPtr = outFilePtr + 'chr' + str(12) + '_1mb_matrix_converg.txt'

swarmCount = 20 # Number of swarms in system
ittCount = 20000 # Maximum itterations before stop
threshold = 0.1 # Error threshold before stoping
randRange = 1.0 # Range of x,y,z starting coords. Random value bewtween -randRange,randRange 

rangeSpace = [] # Max scaling factor. Needs to be optimized for each specific dataset. Use two values [one, two] to multithread through a range of those two values at a interval of 5000


# Arguments for running program
# python3 ParticleChromo3D.py <input_data> <other_parameter>
parser = argparse.ArgumentParser("ParticleChromo3D")
parser.add_argument("infile", help="Matrix of contacts", type=str)
parser.add_argument("-o","--outfile", help="File to output pdb model [Default ./]", type=str, default="./chr.pdb")
parser.add_argument("-rmin","--rangeMin", help="Minimum range for scaling factor. (Input only a min range to do no multithreading) [Default 20000]", type=int, default=20000)
parser.add_argument("-rmax","--rangeMax", help="Minimum range for scaling factor. (If a min and max range are specified, program is multithreaded through intervals of 5000)", type=int)

parser.add_argument("-sc","--swarmCount", help="Number of swarms in system [Default 20]", type=int, default=20)
parser.add_argument("-itt","--ittCount", help="Maximum itterations before stop [Default 20000]", type=int, default=20000)
parser.add_argument("-t","--threshold", help="Error threshold before stoping [Default 0.1]", type=float, default=0.1)
parser.add_argument("-rr","--randRange", help="Range of x,y,z starting coords. Random value bewtween -randRange,randRange [Default 1]", type=float, default=1.0)

args = parser.parse_args()

if args.infile:
	inFilePtr = args.infile
if args.outfile:
	outFilePtr = args.outfile
if args.rangeMin:
	rangeSpace.append(args.rangeMin)
if args.rangeMax:
	rangeSpace.append(args.rangeMax)
if args.swarmCount:
	swarmCount = args.swarmCount
if args.ittCount:
	ittCount = args.ittCount
if args.threshold:
	threshold = args.threshold
if args.randRange:
	randRange = args.randRange

if len(rangeSpace) == 0:
	rangeSpace.append(20000)

if len(rangeSpace) > 2 and (rangeSpace[0] == rangeSpace[1]):
	rangeSpace.pop()

Full_List(rangeSpace, inFilePtr,outFilePtr)