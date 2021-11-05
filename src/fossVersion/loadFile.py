import csv
import numpy as np

class loadFile:
    def __init__(self, file=""):
        ''' Constructor for this class. '''
        self.inFile = file
    
    def matrix2tuple(self, cont):
        ## Scripts convert an input squarematrix to Tuple for mat
        thisLen = len(cont)
        c= []
        for row in range(thisLen):
            for col in  range(row+1,thisLen):
                IF = cont[row][col];       
                c.append([row, col, IF])
        ret = c;
        return ret
    
    def loadFileFunc(self, inFile):
        with open(inFile, 'r') as f:
            reader = csv.reader(f,delimiter='\t')
            cont = []
            for row in reader:
                cont.append(row)
            lenOfRow = len(row)
            if lenOfRow == 3:
                print("this is a sparse matrix")
                stepsize =float(cont[1][1])
            else:
                print("square mat")
                cont = self.matrix2tuple(cont)   
                stepsize = 1
                
        # remove 0 contacts
        tmpCont = []
        for i in range(len(cont)):
            if (float(cont[i][2]) > 0.5):
                tmpCont.append(cont[i][:])
        cont = np.array(tmpCont)
        
        n= len(cont)
        print("Number of points : ", n)
        
        retWithDiag = np.zeros((len(cont), 3))
        
        #TO DO sort
        
        for i in range(len(cont)):
            for j in range(2):
                retWithDiag[i][j] = float(cont[i][j])/stepsize
            retWithDiag[i][2] = float(cont[i][2])
        
        numBins = max(max(retWithDiag[:,0]),max(retWithDiag[:,1]))
        
        print("numBins ",numBins)
        ret = []#np.zeros((int(len(cont)-numBins), 3))

        nextIn = 0
        for i in range(len(cont)-1):
            if retWithDiag[i][0] != retWithDiag[i][1]:
                #ret[nextIn][:] = retWithDiag[i][:]
                ret.append(retWithDiag[i][:])
                nextIn += 1
        
        ret = np.array(ret)
        ret = ret[ret[:,0].argsort()]
        print(ret)
        
        return ret
    

  
