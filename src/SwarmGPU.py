from Swarm import Swarm
import Helper
import copy
import numpy as np
from scipy.spatial import distance
from numba import jit, cuda

class SwarmGPU(Swarm):
    
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
    @jit(target=cuda)
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
    @jit(target=cuda)
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
    @jit(target=cuda)
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