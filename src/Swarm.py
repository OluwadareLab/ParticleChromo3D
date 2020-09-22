import Helper
import copy

import numpy as np

from scipy.spatial import distance

class Swarm:

	id = 0

	def __init__(self, ref, pointCount, randVal=0.5, swarmCount = 10, zeroInd=None, swarmComb = None):

		Swarm.id += 1
		self.id = Swarm.id
		
		self.pc = pointCount
		self.randMax = randVal
		self.randMin = -randVal

		self.gBest = None

		self.ref = ref
		self.zeroInd = zeroInd

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

		self.pos = np.asarray(tempList)
		self.posBest = copy.copy(self.pos)
		self.costBest = np.full((self.pos.shape[0],1), np.inf)

		self.vel = np.zeros((self.pos.shape))
		self.cost = np.full((self.pos.shape[0],1), np.inf)

		self.dist = np.zeros((self.pos.shape[0],self.ref.shape[0]))

		self.locOpCount = np.zeros((self.pos.shape[0],1))

		self.Cost()


	def Rand_Cur(self):
		return np.random.uniform(self.randMin,self.randMax, size=(self.pc, 3))

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

	def Update_Pos(self, itt):
		'''X(t+1) = X(t)+V(t+1)'''

		cutSize = np.random.randint(1, (int)(self.pc-1))

		'''if (self.gBest[1] < 6000):
									thresh = 0.01
								elif (self.gBest[1] < 10000):
									thresh = 0.1
								elif (self.gBest[1] < 25000):
									thresh = 0.5
								else:
									thresh = 1'''

		'''if (itt > 4000):
									thresh = 0.01
								elif(itt > 2000):
									thresh = 0.1
								elif(itt > 500):
									thresh = 0.5'''
		if (itt > 500):
			thresh = (1/itt)*100
		else:
			thresh = 1

		tmpMsk = self.locOpCount > self.Calc_Const(10000,itt,5,15)
		tmpMsk = tmpMsk.reshape(self.locOpCount.shape[0])

		changeInd = np.where(tmpMsk)[0]

		for i in changeInd:
			if(itt < 1000):
				self.pos[i] = self.Rand_Cur()
			else:				
				self.pos[i], mask = self.Rand_shift(self.pos[i], cutSize, thresh)
				self.vel[i][mask] = np.zeros((self.vel.shape[1], self.vel.shape[2]))[mask]
				
		self.vel[tmpMsk] = np.zeros(self.vel.shape)[tmpMsk]
		self.locOpCount[tmpMsk] = -1

		self.pos[~tmpMsk] = self.pos[~tmpMsk]+self.vel[~tmpMsk]

	def Calc_Dist(self):
		for i in range(self.pos.shape[0]):
			if self.zeroInd is None:
				self.dist[i] = distance.pdist(self.pos[i])
			else:
				tmpDist = distance.pdist(self.pos[i])
				tmpDist = np.delete(tmpDist, self.zeroInd)
				self.dist[i] = tmpDist
		
		return self.dist

	# Root Mean Squared Error
	def Cost(self):
		self.Calc_Dist()

		#newCost = np.sqrt( (1/self.pc) * np.sum( (self.dist-ref[:,3])**2, axis=1 ) )
		newCost = np.sum( (self.dist-self.ref[:,3])**2, axis=1 )
		newCost = newCost.reshape((self.pos.shape[0],1))

		self.Update_Cost(newCost)

		return self.cost


	def Calc_Vel(self, ittMax, itt):
		'''V(t+1) = weight*V(t) + conP*ranP(pBest-X(t)) + conG*ranG(gBest-x(t))'''

		ranP = np.random.rand(self.pos.shape[0],self.pos.shape[1],self.pos.shape[2])
		ranG = np.random.rand(self.pos.shape[0],self.pos.shape[1],self.pos.shape[2])

		#weight = self.Calc_Const(ittMax,itt,0.1,0.5)
		#conP = self.Calc_Const(ittMax,itt,0.1,0.5)
		#conG = self.Calc_Const(ittMax,itt,0.1,2)
		weight = 0.5
		conP=0.3
		conG=2.2

		self.vel = (weight*self.vel) + (conP*ranP*(self.posBest-self.pos)) + (conG*ranG*(self.gBest[0]-self.pos))

		return self.vel 

	def Calc_Const(self, ittMax, k, W_MIN, W_MAX):
		if ittMax < k:
			k = ittMax
		return ( ( (W_MAX-W_MIN)*((ittMax-k)/ittMax) ) + W_MIN )