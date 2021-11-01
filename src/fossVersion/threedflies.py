#!/usr/bin/python

import sys
#util
import os
import time

#Math/arrays
import numpy as np
import pandas as pd

from importlib import reload

# Plotting
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d 

#File io
#from Bio.PDB import *
from loadFile import *
from makePDB import *

#PSO
from threeDMaxGrad import *
from myUtils import *
# Import PySwarms
import pyswarms as ps
import numpy.random as random

from importlib import reload
from multiprocessing import Pool

#import pyswarmsDV 

#==========================================================================
# Make directory if it doesn't exist
#--------------------------------------------------------------------------
#if not os.path.exists('my_folder'):
#    os.makedirs('my_folder')

## constants
print ('First argument:',  str(sys.argv[1]))
inFile = str(sys.argv[1])

## load file 
lstCons=loadFile.loadFileFunc(loadFile(), inFile)

lstCons[:,0] = lstCons[:,0].astype('int')
lstCons[:,1] = lstCons[:,1].astype('int')

dataset=lstCons[:,2]
lstCons[:,2] = (dataset - np.min(dataset)) / (np.max(dataset) - np.min(dataset))+0.0001

### config
import yaml
config = yaml.safe_load(open("settings.yml"))
swarm_size = int(config['swarm'][0]['swarm_size'])

al = config['dist_conv'][0]['alpha_low']
ah = config['dist_conv'][0]['alpha_high']
aStepSize = config['dist_conv'][0]['alpha_step_size']
CONVERT_FACTOR_R = np.arange(al, ah, aStepSize) # this is alpha and can be looped

c1 = config['swarm'][0]['c1']
c2 = config['swarm'][0]['c2']
w = config['swarm'][0]['w']
options = {'c1': c1, 'c2': c2, 'w': w}

MAX_ITERATION = config['swarm'][0]['max_iterations']; # maximum number of iterations
topoModel = config['optimization'][0]['model']
p = config['optimization'][1]['extra_params']['p']
k = config['optimization'][1]['extra_params']['k']

threads = int(config['performanceOptions'][0]['numberOfThreads'])
print(threads)

class outputObj:
    def __init__(self, xyzData, outputFile):
        self.xyzData = xyzData
        self.outputFile = outputFile
        
        #defaults
        self.recordName = []
        for i in range(len(xyzData)):
            self.recordName.append('ATOM')
            
def getSpear():
    SUM = 0.0;

    Len = int(n * (n - 1) / 2) 
    Dist = np.zeros(Len);
    WishDist = np.zeros(Len);
    count = 1;
    structure = variables;
    for k in range(len(lstCons)):
        i = int(lstCons[k,0]);    j = int(lstCons[k,1]);    IF = lstCons[k,2];  dist = lstCons[k,3];
        # structure distance   
        x1=structure[i][0];  x2=structure[j][0];
        y1=structure[i][1];  y2=structure[j][1];
        z1=structure[i][2];  z2=structure[j][2];

        str_dist = threeDMaxGrad.calEuclidianDist(threeDMaxGrad(), x1,y1,z1,x2,y2,z2 );
        SUM = SUM + ((str_dist - dist)**2);

        # calculate spearman_correlation and Pearson correlation
        if (i != j and IF > 0 and count < Len  ):   
            Dist[count] = str_dist;
            WishDist[count]= dist;
            count = count + 1

    SUM = SUM / len(lstCons);    
    rmse = np.sqrt(SUM);

    # let's convert to a dataframe
    df = pd.DataFrame({'Dist': Dist, 'wishDist': WishDist})
    pearsonCoeff = df.corr(method = 'pearson')
    spearmanCoeff = df.corr(method = 'spearman')
    print("this spear  : ", spearmanCoeff.iloc[0]['wishDist']) 
    
def opt_func(X):
    n_particles = X.shape[0]  # number of particles
    
    d = 0;
    #(n)
    #print(max(X[0,:]))
    #print(min(X[0,:]))

    #structure = variables #this is xyz
    allDists = np.zeros(n_particles)
    
    cost = 0
    for particle in range(n_particles):
        for i in range(len(structure)):
            structure[i][0]=X[particle][i]
            structure[i][1]=X[particle][i+len(structure)]
            structure[i][2]=X[particle][i+len(structure)*2]
            
        
        #loop through existing data
        for k in range(len(lstCons)):
            #IF = lstCons[k,2];  
            if (lstCons[k,2] <= 0) :
                continue
            i = int(lstCons[k,0]);    
            j = int(lstCons[k,1]);    
            
            dist = lstCons[k,3];
            
            
                
            point1 = np.array((structure[i][0], structure[i][1], structure[i][2]))
            point2 = np.array((structure[j][0], structure[j][1], structure[j][2]))

            str_dist = np.linalg.norm(point1 - point2)     

            # objective function  
            d += ((str_dist - dist)**2)
            
            #cost
            #cost += -(n/2) - (n*np.log(np.sqrt(d/n)));
            
        
        allDists[particle] = d/len(structure)
    
    return allDists
    
def opt_func_alt(X):
    n_particles = X.shape[0]  # number of particles
    
    d = 0;
    #print(X)


    structure = variables #this is xyz
    allDists = np.zeros(n_particles)
    
    cost= 0
    
    for particle in range(n_particles):
        for i in range(len(structure)):
            structure[i][0]=structure[i][0]+X[particle][i]
            structure[i][1]=structure[i][1]+ X[particle][i+len(structure)]
            structure[i][2]=structure[i][2]+X[particle][i+len(structure)*2]
            
        
        #loop through existing data
        for k in range(len(lstCons)):
            #IF = lstCons[k,2];  
            if (lstCons[k,2] <= 0) :
                continue
            i = int(lstCons[k,0]);    
            j = int(lstCons[k,1]);    
            
            dist = lstCons[k,3];
            
            
                
            point1 = np.array((structure[i][0], structure[i][1], structure[i][2]))
            point2 = np.array((structure[j][0], structure[j][1], structure[j][2]))

            str_dist = np.linalg.norm(point1 - point2)     

            # objective function  
            d += ((str_dist - dist)**2)
            
            #cost
            #cost += -(n/2) - (n*np.log(np.sqrt(d/n)));
            
        
        allDists[particle] = d/len(structure)
    
    return allDists
    
def opt_func_oneGuy(X):
    n_particles = X.shape[0]  # number of particles
    
    d = 0;
    #print(X)


    structure = variables #this is xyz
    allDists = np.zeros(n_particles)
    
    cost= 0
    
    for particle in range(n_particles):

        structure[thisXYZ][0]=structure[thisXYZ][0] + X[particle][0]
        structure[thisXYZ][1]=structure[thisXYZ][0] + X[particle][1]
        structure[thisXYZ][2]=structure[thisXYZ][0] + X[particle][2]
            
        
        #loop through existing data
        for k in range(len(lstCons)):
            #IF = lstCons[k,2];  
            if (lstCons[k,2] <= 0) :
                continue
                
            
            i = int(lstCons[k,0]);    
            j = int(lstCons[k,1]);    
            if (i != thisXYZ):
                continue

            dist = lstCons[k,3];

            point1 = np.array((structure[i][0], structure[i][1], structure[i][2]))
            point2 = np.array((structure[j][0], structure[j][1], structure[j][2]))

            str_dist = np.linalg.norm(point1 - point2)     

            # objective function  
            d += ((str_dist - dist)**2)
            
            #cost
            #cost += -(n/2) - (n*np.log(np.sqrt(d/n)));
            
        
        allDists[particle] = d/len(structure)
    
    return allDists


AVG_DIST = 10  # an arbitrary distance

n=int(max(max(lstCons[:,0]),max(lstCons[:,1])))+1
print("n = ",n)

bestAlpha= CONVERT_FACTOR_R[0]
lstConsReset = lstCons
for CONVERT_FACTOR in CONVERT_FACTOR_R : 
    lstCons = lstConsReset
    ## Find the average IF
    avgIF = 0.0
    for i in range(len(lstCons)):
        avgIF = avgIF + float(lstCons[i][2])
    avgIF = avgIF/len(lstCons)
    
    maxIF = 0.0
    ## scale average distance to AVG_DIST
    avgDist = AVG_DIST;
    avgAdjIF = 0.0;
    avgAdjCount = 0;
    totalIF = 0;

    for i in range(len(lstCons)):
        x = lstCons[i][0]
        y = lstCons[i][1]
        IF = lstCons[i][2]
        lstCons[i][2] = IF / avgIF  # normallize IF by avgIF
        IF = lstCons[i][2]
        dist = 1/(IF**CONVERT_FACTOR)
        avgDist = avgDist + dist

        totalIF = totalIF +  IF

        if ( IF > maxIF):
            maxIF = lstCons[i][2]
        # Find the adjacent position IF
        if ( abs(x-y)==1):
            avgAdjCount= avgAdjCount+1
            avgAdjIF =  avgAdjIF + IF

    avgDist = avgDist/len(lstCons)
    avgAdjIF = avgAdjIF/avgAdjCount

    maxIF = min(avgAdjIF, maxIF)

    ## TODO Add adjacent if none exist
    
    print('TODO Added missing adjacent constraint...')

    print('Number of constraints: = ', n)
    maxD = 0
    distResultsList= []
    for i in range(len(lstCons)):
        IF = lstCons[i,2];
        dist = AVG_DIST/ ((IF**CONVERT_FACTOR)* avgDist)
        distResultsList.append(dist)
        if (dist > maxD):
            maxD = dist;
    #print("rand dists are : ", sum(distResultsList)/len(distResultsList))
    result = np.hstack((lstCons, np.atleast_2d(distResultsList).T))
    print(lstCons)
    lstCons = result
    print('Max distance is: = \n', maxD); 
    
    ## Optimization
    ## Initialize Structure
    #=========================================================================

    thisStr =  [];
    R = [0,12.5];
    for i in range(n):
        xyz = np.array([random.random(),random.random(),random.random()]) * (R[1]-R[0]) + min(R) 
        thisStr.append(xyz)
    
    ## Variables declaration
    #=========================================================================
    thislen = len(lstCons)

    variables = thisStr;
    oldobj = 0;
    getSpear()
    ######### PSO
    variables = np.array(variables)

    dist = 1/(IF**CONVERT_FACTOR)
    dim = n*3       # Dimension of X
    epsilon = 1
    
    
    max_bound = maxD/2 * np.ones(dim)
    min_bound = - max_bound
    bounds = (min_bound, max_bound)
    
    
    lr = 0.5
    clamp = (-lr , lr)
    
    if (topoModel == "pyramid"):
        from pyswarms.backend.topology import Pyramid
        my_topology = Pyramid(static=False)
        optimizer = ps.single.GeneralOptimizerPSO(n_particles=swarm_size,
                                        dimensions=dim,
                                        options=options , velocity_clamp = clamp, topology=my_topology)
    elif (topoModel == "random"):
        from pyswarms.backend.topology import Random as PSOTopoRandom
        my_topology = PSOTopoRandom()
        options = {'c1': c1, 'c2': c2, 'w': w, 'k': k, 'p':p}
        optimizer = ps.single.GeneralOptimizerPSO(n_particles=swarm_size,
                                        dimensions=dim,
                                        options=options , velocity_clamp = clamp, topology=my_topology)
    elif (topoModel == "star"):
        from pyswarms.backend.topology import Star
        my_topology = Star(static=False)
        optimizer = ps.single.GeneralOptimizerPSO(n_particles=swarm_size,
                                        dimensions=dim,
                                        options=options , velocity_clamp = clamp, topology=my_topology)
    elif (topoModel == "local"):
        options = {'c1': c1, 'c2': c2, 'w': w, 'k': k, 'p': p}
        optimizer = ps.single.LocalBestPSO(n_particles=swarm_size,
                                        dimensions=dim,
                                        options=options , velocity_clamp = clamp)
    else:
        print("Chose global best with ", swarm_size, " dims : ", dim, "options : ", options)
        print()
        optimizer = ps.single.GlobalBestPSO(n_particles=swarm_size,
                                        dimensions=dim,
                                        options=options , velocity_clamp = clamp)
    
    
    # Perform optimization
    print("here")
    structure = variables #this is xyz
    
    
    
    cost, bestPXYZ = optimizer.optimize(opt_func, iters=MAX_ITERATION, n_processes=threads , verbose=bool(config['performanceOptions'][1]['verbose'])==1)   
    for i in range(len(variables)):
        variables[i][0] = bestPXYZ[i]
        variables[i][1] = bestPXYZ[i+len(variables)]
        variables[i][2] = bestPXYZ[i+len(variables)*2]
    #print("initRunCost : ", cost)
    getSpear()
    '''
    lr = 1
    lastCost = cost
    options = {'c1': 0.9, 'c2':0.9, 'w':0.9, 'k': 3, 'p': 2}#Local Best
    noImprovment = 0
    for i in range(n):
        optimizer = ps.single.LocalBestPSO(n_particles=swarm_size,
                                    dimensions=3,
                                    options=options )
        thisXYZ =  int(random.random()*n)
        
        thisCost, bestPXYZ = optimizer.optimize(opt_func_oneGuy, iters=50, n_processes=10, verbose=False)   
        
        if lastCost-thisCost > 0:#Less then 1 % improvment
            #for i in range(len(variables)):
            variables[thisXYZ][0] += lr*bestPXYZ[0]
            variables[thisXYZ][1] += lr*bestPXYZ[1]
            variables[thisXYZ][2] += lr*bestPXYZ[2]
            lastCost=thisCost
        else: 
            print("in else")
            noImprovment += 1
            if(noImprovment == n):
                lr = lr / 10'''
            
    #========================================================================
    # scoring using spearman correlation, pearson correlation and  RMSD     
    #------------------------------------------------------------------------
    # calculate rmse    
    SUM = 0.0;

    Len = int(n * (n - 1) / 2) 
    Dist = np.zeros(Len);
    WishDist = np.zeros(Len);
    count = 1;
    structure = variables;
    for k in range(len(lstCons)):
        i = int(lstCons[k,0]);    j = int(lstCons[k,1]);    IF = lstCons[k,2];  dist = lstCons[k,3];
        # structure distance   
        x1=structure[i][0];  x2=structure[j][0];
        y1=structure[i][1];  y2=structure[j][1];
        z1=structure[i][2];  z2=structure[j][2];

        str_dist = threeDMaxGrad.calEuclidianDist(threeDMaxGrad(), x1,y1,z1,x2,y2,z2 );
        SUM = SUM + ((str_dist - dist)**2);

        # calculate spearman_correlation and Pearson correlation
        if (i != j and IF > 0 and count<Len  ):   
            Dist[count] = str_dist;
            WishDist[count]= dist;
            count = count + 1

    SUM = SUM / len(lstCons);    
    rmse = np.sqrt(SUM);

    # let's convert to a dataframe
    df = pd.DataFrame({'Dist': Dist, 'wishDist': WishDist})
    pearsonCoeff = df.corr(method = 'pearson')
    spearmanCoeff = df.corr(method = 'spearman')
   
    if CONVERT_FACTOR == CONVERT_FACTOR_R[0]:#first run
        bestSpearmanRHO = spearmanCoeff.iloc[0]['wishDist']
        bestPearsonRHO = pearsonCoeff.iloc[0]['wishDist']
        bestMat =  variables
    elif bestSpearmanRHO < spearmanCoeff.iloc[0]['wishDist']:
        bestSpearmanRHO = spearmanCoeff.iloc[0]['wishDist']
        bestPearsonRHO = pearsonCoeff.iloc[0]['wishDist']
        bestMat =  variables
        bestAlpha = CONVERT_FACTOR
#========================================================================
# create pdb
#------------------------------------------------------------------------
#Increase structure Size

xyz4pdb = myUtils.convert2xyz(myUtils(), n, bestMat) 
scale=100/np.amax(xyz4pdb)
xyz4pdb = xyz4pdb* scale
#output pdb file
outputData= np.round_(xyz4pdb ,3 )

outputFile = 'yolo.pdb' #output directory.


if os.path.exists(outputFile):
    os.remove(outputFile)

output = outputObj(outputData,outputFile)

pdbMaker = makePDB()
makePDB.mat2pdb(pdbMaker, output) # Converts the mat XYZ coordinates to PDB format.

print("Input file: ", inFile)
print("Convert factor:: ",bestAlpha)
print("AVG RMSE  : ", rmse)    
print("AVG Spearman correlation Dist vs. Reconstructed Dist  : ", bestSpearmanRHO) 
print("AVG Pearson correlation Dist vs. Reconstructed Dist  : ", bestPearsonRHO) 


FILE = open("this.log", 'w');
inString = ("Input file: " + str(inFile) + "\n"+"Convert factor:: " + str(bestAlpha)+ "\n" +
    "AVG RMSE  : "+ str(rmse)+ "\n"+"AVG Spearman correlation Dist vs. Reconstructed Dist  : "+ str(bestSpearmanRHO)+ "\n"
    +"AVG Pearson correlation Dist vs. Reconstructed Dist  : "+ str(bestPearsonRHO)) + "\n"
    
FILE.write(inString);
FILE.close();