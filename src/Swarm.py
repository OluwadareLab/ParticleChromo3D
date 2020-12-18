import Helper
import copy

import numpy as np

from scipy.spatial import distance


''' Swarm Cluster for PSO calculations on HiC data to find optimal distances for xyz points of beads
        This creates a cluster of swarms, each swarm is represented by a list of xyz coordinates, 
        a particle of each swarm is the individual coordinate (x,y, or z) of each bead. A different velocity is applied
        to each individual coordinate (x, y, or z).

    Swarms are in a cluster format as this allows python to do vector math on the swarm array which is MANY times more efficient
    than doing individual objects of swarms
    
    Developer: 
    David Vadnais dvadnais@uccs.edu
    DMichael Middleton mmiddlet@uccs.edu 12/4/2020

'''

class Swarm:

    id = 0

    # ref: distance matrix from HiC data
    # pointCount: Number of beads
    # randVal: random value to calculate initial x,y,z from
    # swarmCount: number of swarms in cluster
    # zeroInd: Used to delete in distance calculation
    # swarmComb: used to combine multiple swarms if doing multiple passes
    def __init__(self, ref, pointCount, randVal=0.5, swarmCount = 10, zeroInd=None, swarmComb = None):

        Swarm.id += 1
        self.id = Swarm.id # ID of swarm if using multiprocession with multiple swarm clusters

        self.pc = pointCount # Pointcount

        # Random Values to get initial xyz coordinates from
        self.randMax = randVal
        self.randMin = -randVal

        self.gBest = None # Global best position

        self.ref = ref # Reference distance matrix
        self.zeroInd = zeroInd

        # Creates a list of swarms where each swarm is a list of xyz coordinates of size pointcount
        tempList = []
        if swarmComb is None:
            for i in range(swarmCount):
                tempList.append(self.Rand_Cur())			
        else:
            for swarm in swarmComb:
                tempList.append(swarm.gBest[0])

            for i in range(swarmCount-len(swarmComb)):
                cutSize = np.random.randint(1, (int)(self.pc/2))
                randCopy = np.random.choice(swarmComb).gBest[0]
                tempPos, mask = self.Rand_shift(randCopy, cutSize, 0.1)
                tempList.append(tempPos)


        self.pos = np.asarray(tempList) # Turns the list into a 3D matrix for vector operations
        self.posBest = copy.copy(self.pos) # Best position of each individual swarm
        self.costBest = np.full((self.pos.shape[0],1), np.inf) # Best cost of each individual swarm

        self.vel = np.zeros((self.pos.shape)) # Velocity of each particle in each swarm
        self.cost = np.full((self.pos.shape[0],1), np.inf) # Cost of each swarm

        self.dist = np.zeros((self.pos.shape[0],self.ref.shape[0])) # Distance of each swarm

        self.locOpCount = np.zeros((self.pos.shape[0],1)) # Used to check if a swarm is not changing after a certain amount of itterations

        self.Cost() # Gets first cost calculations

    # Gets a array of xyz coordinates of size pointcount,3 between a minimum and maximum value
    def Rand_Cur(self):
        return np.random.uniform(self.randMin,self.randMax, size=(self.pc, 3))

    # Performs a random shift of all the xyz positions in a passed in position array
    # copyPos: position array to shift
    # cutSize: Number of values to shift
    # Threshold: size of shift
    def Rand_shift(self, copyPos, cutSize, threshold=0.1):

        temp = copy.copy(copyPos)

        # random boolean mask for which values will be changed
        falseVal = np.zeros(cutSize).astype(np.bool)
        trueVal = np.ones(temp.shape[0]-cutSize).astype(np.bool)

        randBool = np.append(falseVal, trueVal)
        np.random.shuffle(randBool)

        mask =  np.tile(randBool.transpose(), (temp.shape[1], 1)).transpose()

        # random matrix the same shape of your data
        r = np.random.uniform(-threshold, threshold, size=temp.shape)

        temp[mask] += r[mask]

        return temp, mask

    # Updates the cost of each swarm in cluster
    def Update_Cost(self, newCost):
        tmpMsk = newCost > self.cost

        self.locOpCount[tmpMsk] += np.ones(self.locOpCount.shape)[tmpMsk]

        self.cost = newCost

        tmpMsk = self.cost < self.costBest
        tmpMsk = tmpMsk.reshape(tmpMsk.shape[0])

        self.posBest[tmpMsk] = copy.copy(self.pos)[tmpMsk]
        self.costBest[tmpMsk] = self.cost[tmpMsk]

        currentBest = np.argmin(self.cost);

        if (self.gBest is None) or (self.cost[currentBest] < self.gBest[1]):
            self.gBest = (copy.copy(self.pos[currentBest]),copy.copy(self.cost[currentBest][0]), copy.copy(self.dist[currentBest]))

    # Updates the position of each swarm
    # itt: the current itteration of program
    def Update_Pos(self, itt):
        '''X(t+1) = X(t)+V(t+1)'''

        # Gets a random cutize
        cutSize = np.random.randint(1, (int)(self.pc-1))

        # Adjust threshold depending on itteration to gradually adjust
        if (itt > 500):
            thresh = (1/itt)*100
        else:
            thresh = 1

        # Gets a boolean array of values that checks if a certin number of particles has changed in the last itteration
        # This certain number is calculated in Calc_Const and is based on the number of itterations
        # The boolean array that is calculated is then used to shift any swarms that have been calculated to not been changing
        tmpMsk = self.locOpCount > self.Calc_Const(10000,itt,5,15)
        tmpMsk = tmpMsk.reshape(self.locOpCount.shape[0])
        changeInd = np.where(tmpMsk)[0]

        # For each swarm that is being changed
        for i in changeInd:
            # If itterations are over a certain amount a full new swarm is calculated
            if(itt < 1000):
                self.pos[i] = self.Rand_Cur()
            else:	
                # Otherwise, shift only the cutsize amount of particles			
                self.pos[i], mask = self.Rand_shift(self.pos[i], cutSize, thresh)
                self.vel[i][mask] = np.zeros((self.vel.shape[1], self.vel.shape[2]))[mask]

        # Changes all shifted particles to have a velocity of zero
        self.vel[tmpMsk] = np.zeros(self.vel.shape)[tmpMsk]
        # Resets local optima counter for changed swarms
        self.locOpCount[tmpMsk] = -1

        # Adjusts new positions with velocitys
        self.pos[~tmpMsk] = self.pos[~tmpMsk]+self.vel[~tmpMsk]

    # Calculates the euclidean distance for each swarm
    def Calc_Dist(self):
        for i in range(self.pos.shape[0]):
            if self.zeroInd is None:
                self.dist[i] = distance.pdist(self.pos[i])
            else:
                tmpDist = distance.pdist(self.pos[i])
                tmpDist = np.delete(tmpDist, self.zeroInd)
                self.dist[i] = tmpDist

        return self.dist
    
    # Calculates the cost of each swarm
    # Sum Squared Error
    def Cost(self):
        self.Calc_Dist()

        #newCost = np.sqrt( (1/self.pc) * np.sum( (self.dist-ref[:,3])**2, axis=1 ) )
        newCost = np.sum( (self.dist-self.ref[:,3])**2, axis=1 )#SSE
        #newCost = np.sqrt(np.sum( (self.dist-self.ref[:,3])**2, axis=1 ))#RMSE
        #delta = 1.0
        #y = self.ref[:,3]
        #yHat = self.dist
        #newCost = np.sum(np.where(np.abs(y-yHat) < delta,.5*(y-yHat)**2 , delta*(np.abs(y-yHat)-0.5*delta)),axis=1)
        #print("ME:" ,newCost)
        #print("notME",dnewCost)
        
        newCost = newCost.reshape((self.pos.shape[0],1))

        self.Update_Cost(newCost)

        return self.cost

    # Calculates the velocitys of each particle
    def Calc_Vel(self, ittMax, itt):
        '''V(t+1) = weight*V(t) + conP*ranP(pBest-X(t)) + conG*ranG(gBest-x(t))'''

        ranP = np.random.rand(self.pos.shape[0],self.pos.shape[1],self.pos.shape[2])
        ranG = np.random.rand(self.pos.shape[0],self.pos.shape[1],self.pos.shape[2])

        # Way to calculate constants, not used anymore but for future reference
        #weight = self.Calc_Const(ittMax,itt,0.1,0.5)
        #conP = self.Calc_Const(ittMax,itt,0.1,0.5)
        #conG = self.Calc_Const(ittMax,itt,0.1,2)

        # Constants found to work the best
        weight = 0.5
        conP=0.3
        conG=2.2

        self.vel = (weight*self.vel) + (conP*ranP*(self.posBest-self.pos)) + (conG*ranG*(self.gBest[0]-self.pos))

        return self.vel 

    # Function that calculates a constant value
    # ittMax: maximum number of itterawtions
    # k: current itterations
    # W_MIN: constant min
    # W_MAX: constant max
    def Calc_Const(self, ittMax, k, W_MIN, W_MAX):
        if ittMax < k:
            k = ittMax
        return ( ( (W_MAX-W_MIN)*((ittMax-k)/ittMax) ) + W_MIN )