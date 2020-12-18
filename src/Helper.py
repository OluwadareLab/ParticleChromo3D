import Swarm

from Bio import PDB
from scipy.spatial import procrustes

from collections.abc import Iterable

import numpy as np
import copy

# Reads in a matrix of contacts
def Read_Matrix_To_List(filePtr):
    contact = np.genfromtxt(filePtr, delimiter=' ')
    #delete zero contacts
    contact = contact[~np.all(contact == 0, axis=0)]
    idx = np.argwhere(np.all(contact[..., :] == 0, axis=0))
    contact = np.delete(contact, idx, axis=1)
    
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

# Reads in a list format for contacts
def Read_Data_List(filePtr):
    return np.loadtxt(filePtr), None

# Reads the data and puts it into a matrix of distances
def Read_Data(filePtr, convFactor=None):
    constraint, zeroInd = Read_Matrix_To_List(filePtr)
    #constraint, zeroInd = Read_Data_List(filePtr)

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

    mean = np.mean(constraint[:,2])

    #dist = np.zeros(constraint.shape[0])
    #if maxScale is not None:
        #constraint[:,2] = Scale_Arr(constraint[:,2],1000,20000)
        #constrAvg = avgCalc(constraint, convFactor)
    #dist = ( 1.0 / (constraint[:,2]**convFactor) )

    #constraint = np.insert(constraint,3,dist,axis=1)

    return constraint, pointMap, zeroInd

# Calcualtes the average number of contacts
def avgCalc(constraint, convFactor):
    avgIf = constraint.mean(axis=0)[2]
    avgDist = 0

    constrAvg = copy.copy(constraint[:,2])

    for i in range(constraint.shape[0]):
        constrAvg[i] = constrAvg[i]/avgIf

    #avgDist = (constrAvg**convFactor).mean()

    return constrAvg#, avgDist

# Writes xyz list to a PDB
def Write_Output(filePtr, xyz):
    xyz = Scale_Arr(xyz)
    WritePDB(xyz, filePtr+".pdb")

# Scales a range of values
def Scale_Arr(xyz, minVal=-10, maxVal=10):
    oldRange = xyz.max() - xyz.min()
    newRange = (maxVal - minVal)  
    return (((xyz - xyz.min()) * newRange) / oldRange) + minVal

# Writes a list to a file
def Write_List(listToWrite, filePtr):
    with open(filePtr, 'w') as f:
        for item in listToWrite:
            if isinstance(item, Iterable):
                for subItem in item:
                    f.write("%s " % subItem)
                f.write("\n")
            else:
                f.write("%s\n" % item)

# Performs procrustes on a pdb file
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
            #j = 1
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

def Write_Log(outfile, inFile, bestAlpha, rmse, bestSpearmanRHO,bestPearsonRHO):
    FILE = open(outfile, 'w');
    inString = ("Input file: " + str(inFile) + "\n"+"Convert factor:: " + str(bestAlpha)+ "\n" +
        "Best cost  : "+ str(rmse)+ "\n"+"Best Spearman correlation Dist vs. Reconstructed Dist  : "+ str(bestSpearmanRHO)+ "\n" + "Best Pearson correlation Dist vs. Reconstructed Dist  : "+ str(bestPearsonRHO)) + "\n"

    FILE.write(inString);
    FILE.close();
    
    