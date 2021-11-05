import numpy as np
class myUtils:
    def __init__(self):
        ''' Constructor for this class. '''
        
        
    def convert2xyz(self, n,variables):
        ## convert to xyz
        xyzTmp  = np.zeros((n,3))
        for i in range(n) :
            xyzTmp[i-1,0] = variables[i-1][0]
            xyzTmp[i-1,1] = variables[i-1][1]
            xyzTmp[i-1,2] = variables[i-1][2]
        return xyzTmp
    
    def isconvergence(self, change, cost, NEAR_ZERO):
        ## check if the size of derivatives/gradient is close to zero
        converge = 0;
        SUM = 0;
        sq_change = change**2;

        SUM  = np.sum(sq_change);
        SUM  = np.sqrt(SUM);

        if (SUM < (NEAR_ZERO * abs(cost))):
            converge = 1;
        return converge