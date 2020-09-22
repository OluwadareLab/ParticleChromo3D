import Swarm

from Bio import PDB
from scipy.spatial import procrustes

from collections.abc import Iterable

import numpy as np
import copy

def Read_Matrix_To_List(filePtr):
	contact = np.genfromtxt(filePtr, delimiter=' ')

	print(contact.shape)
	contact = contact[~np.all(contact == 0, axis=0)]
	idx = np.argwhere(np.all(contact[..., :] == 0, axis=0))
	contact = np.delete(contact, idx, axis=1)
	print(contact.shape)

	count = 0

	zeroInd = []
	stopDupe = {}
	contactList = []
	for i in range(contact.shape[0]):
		row =  contact[i]
		for j in range(row.shape[0]):
			if (i != j) and ((i,j) not in stopDupe) and ((j,i) not in stopDupe):
				stopDupe[(i,j)] = 1
				stopDupe[(j,i)] = 1
				if (row[j] > 0):
					contactList.append( (i,j, row[j]) )
				else:
					zeroInd.append(count)
				count += 1

	if len(zeroInd) <= 0:
		zeroInd=None

	return np.array(contactList), zeroInd

def Read_Data_List(filePtr):
	return np.loadtxt(filePtr), None

def Read_Data(filePtr, maxScale, convFactor=None):
	constraint, zeroInd = Read_Matrix_To_List(filePtr)

	pointList = set()
	for i in range(constraint.shape[0]):
		pointList.add(constraint[i][0])
		pointList.add(constraint[i][1])

	pointMap = {}
	for i in range(len(pointList)):
		pointMap[pointList.pop()] = i

	for i in range(constraint.shape[0]):
		constraint[i][0] = pointMap[constraint[i][0]]
		constraint[i][1] = pointMap[constraint[i][1]]

	#constraint[:,2] = constraint[:,2]*50
	mean = np.mean(constraint[:,2])
	print(mean)
	#maxScale = (20000000*(mean**-0.783))

	#maxScale = int((156.05*len(pointMap)) + 3973.3)
	#maxScale = int((117.39*len(pointMap)) + 5323.5)
	#maxScale = int((112.26*len(pointMap)) + 5954.9)

	#print(maxScale)
	#maxScale = int(np.ceil(maxScale / 5000) * 5000)
	#maxScale = int(np.floor(maxScale / 5000) * 5000)
	#maxScale = int(np.round(maxScale / 5000) * 5000)

	dist = np.zeros(constraint.shape[0])
	if maxScale is not None:
		constraint[:,2] = Scale_Arr(constraint[:,2],1000,maxScale)
		constrAvg, avgDist = avgCalc(constraint, convFactor)
		dist = ( 10 / ((constrAvg**convFactor)*avgDist) )
	else:
		saveLow = float('inf')
		stepVals = range(2000, 100000, 2000)
		for i in stepVals:
			tmpConstraint = copy.copy(constraint)
			tmpConstraint[:,2] = Scale_Arr(tmpConstraint[:,2],1000, i)

			constrAvg, avgDist = avgCalc(tmpConstraint, convFactor)
			tmpDist = ( 10 / ((constrAvg**convFactor)*avgDist) )

			tmpConstraint = np.insert(tmpConstraint,3,tmpDist,axis=1)
			swarm = Swarm.Swarm(tmpConstraint, len(pointMap), randVal=1, swarmCount=5, zeroInd=zeroInd)

			if np.abs(swarm.gBest[1]) < saveLow:
				dist = copy.copy(tmpDist)
				saveLow = np.abs(swarm.gBest[1])
				saveRange = i
				print(saveRange)

	'''	if len(pointMap) > 600:
			constraint[:,2] = Scale_Arr(constraint[:,2],1000,100000)
		elif len(pointMap) > 260:
			constraint[:,2] = Scale_Arr(constraint[:,2],1000,50000)
		elif len(pointMap) > 200:
			constraint[:,2] = Scale_Arr(constraint[:,2],1000,40000)
		elif len(pointMap) > 150:
			constraint[:,2] = Scale_Arr(constraint[:,2],1000,30000)
		elif len(pointMap) > 50:
			constraint[:,2] = Scale_Arr(constraint[:,2],1000,20000)
		else:
			constraint[:,2] = Scale_Arr(constraint[:,2],1000,10000)'''

	constraint = np.insert(constraint,3,dist,axis=1)

	print(constraint.shape)

	#dist = 1 / (constraint[:,2]**convFactor)

	return constraint, pointMap, zeroInd

def avgCalc(constraint, convFactor):
	avgIf = constraint.mean(axis=0)[2]
	avgDist = 0

	constrAvg = copy.copy(constraint[:,2])

	for i in range(constraint.shape[0]):
		constrAvg[i] = constrAvg[i]/avgIf

	avgDist = (constrAvg**convFactor).mean()

	return constrAvg, avgDist

def Write_Output(filePtr, xyz):
	xyz = Scale_Arr(xyz)
	WritePDB(xyz, filePtr+'.pdb')

def Scale_Arr(xyz, minVal=-10, maxVal=10):
	oldRange = xyz.max() - xyz.min()
	newRange = (maxVal - minVal)  
	return (((xyz - xyz.min()) * newRange) / oldRange) + minVal


def Write_List(listToWrite, filePtr):
	with open(filePtr, 'w') as f:
		for item in listToWrite:
			if isinstance(item, Iterable):
				for subItem in item:
					f.write("%s " % subItem)
				f.write("\n")
			else:
				f.write("%s\n" % item)

def Proc_PDB(inputPtr, comparePtr):
	xyzIn = Read_PDB(inputPtr)
	xyzComp = Read_PDB(comparePtr)
	standXYZComp, newXYZ, dis = procrustes(xyzComp,xyzIn)
	print(xyzIn)
	print(xyzComp)
	print(newXYZ)
	return standXYZComp,newXYZ

'''SOURCE: https://github.com/mbglab/EVR'''
def WritePDB(positions, pdb_file, ctype = "0"):
	'''Save the result as a .pdb file'''
	o_file = open(pdb_file, "w")
	o_file.write("\n")

	col1 = "ATOM"
	col3 = "CA MET"
	col8 = "0.20 10.00"

	bin_num = len(positions)

	for i in range(1, bin_num+1):
		col2 = str(i)
		col4 = "B"+col2
		col5 = "%.3f" % positions[i-1][0]
		col6 = "%.3f" % positions[i-1][1]
		col7 = "%.3f" % positions[i-1][2]
		col2 = " "*(5 - len(col2)) + col2
		col4 = col4 + " " * (6 - len(col4))
		col5 = " " * (8 - len(col5)) + col5
		col6 = " " * (8 - len(col6)) + col6
		col7 = " " * (8 - len(col7)) + col7

		col = (col1, col2, col3, col4, col5, col6, col7,col8)
		line = "%s  %s   %s %s   %s%s%s  %s\n" % col
		o_file.write(line)
	col1 = "CONECT"
	for i in range(1, bin_num+1):
		col2 = str(i)
		j = i + 1
		if j > bin_num:
			if ctype == "1":
				continue
			j = 1
		col3 = str(j)

		col2 = " " * (5 - len(col2)) + col2
		col3 = " " * (5 - len(col3)) + col3

		line = "%s%s%s\n" % (col1, col2, col3)
		o_file.write(line)

	o_file.write("END")
	o_file.close()

'''SOURCE: https://stackoverflow.com/questions/47825542/how-to-separately-get-the-x-y-or-z-coordinates-from-a-pdb-file'''
def Read_PDB(filePtr):
	parser = PDB.PDBParser()
	io = PDB.PDBIO()
	struct = parser.get_structure('temp',filePtr)

	xyz = []
	for model in struct:
	    for chain in model:
	        for residue in chain:
	            for atom in residue:
	                x,y,z = atom.get_coord()
	                xyz.append(np.array([x,y,z]))

	return np.asarray(xyz)	 
