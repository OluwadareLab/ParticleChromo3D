import numpy as np

class threeDMaxGrad:
    def __init__(self):
        ''' Constructor for this class. '''
        
    ## Gradient Calculator
    def gradientCalculator(self, thislen, variables, lstCons, n, maxIF, dist):
        # Initialize
        cost = 0
        change = np.zeros((thislen,3))
        val = 0;
        structure = variables;

        # Calculate the chain rule derivative for the gradient calculation
        tmper = self.gradChainCalc(lstCons,structure,n)
        dl_dw = tmper[0]
        dw_dv = tmper[1]

        #loop through existing data
        for k in range(len(lstCons)):
            i = int(lstCons[k,0]);    j = int(lstCons[k,1]);    IF = lstCons[k,2];  dist = lstCons[k,3];
            if (IF <= 0) :
                continue


            if (abs(i - j) == 1 ) :
                IF = 1.0 * maxIF

            x1=structure[i][0];  x2=structure[j][0];
            y1=structure[i][1];  y2=structure[j][1];
            z1=structure[i][2];  z2=structure[j][2];

            str_dist = self.calEuclidianDist(x1,y1,z1,x2,y2,z2 );     

            z = str_dist - dist;

            # the remaining part of dv_d(x,y,z)

            tmp =   dl_dw * dw_dv * 2 * (z/str_dist);

            # objective function       
            val = val  + (z**2);        


            change[i,0] = change[i,0] + (tmp * (structure[i][0] - structure[j][0]));
            change[i,1] = change[i,1] + (tmp * (structure[i][1] - structure[j][1]));
            change[i,2] = change[i,2] + (tmp * (structure[i][2] - structure[j][2]));

            change[j,0] = change[j,0] + (tmp * (structure[j][0] - structure[i][0]));
            change[j,1] = change[j,1] + (tmp * (structure[j][1] - structure[i][1]));
            change[j,2] = change[j,2] + (tmp * (structure[j][2] - structure[i][2]));



        # cost
        cost = -(n/2) - (n*np.log(np.sqrt(val/n)));

        return change , cost
    
    def gradChainCalc(self, lstCons, structure, n):
        v = 0;
        for k in range(len(lstCons)):
            i = int(lstCons[k,0]);    j = int(lstCons[k,1]);    IF = lstCons[k,2];  dist = lstCons[k,3];
            if (IF <= 0):
                continue;
               # structure distance   
            #print(i,"  ",j)
            x1=structure[i][0];  x2=structure[j][0];
            y1=structure[i][1];  y2=structure[j][1];
            z1=structure[i][2];  z2=structure[j][2];
            str_dist = self.calEuclidianDist(x1,y1,z1,x2,y2,z2 );      
            # dist = IF distance 
            z = ((str_dist-dist)**2);
            v = v + z;



        ## Calculate w
        #------------------------------------------------------------------------
        w = np.sqrt(v/n);
        #=========================
        # Find dl_dw
        #=========================
        dl_dw = -n/w;
        #=========================
        #Find dw_dv
        #=========================
        dw_dv = 1/(2*np.sqrt(n*v))

        return dl_dw, dw_dv

    def calEuclidianDist(self, x1,y1,z1,x2,y2,z2):
        # calEuclidianDist = square euclidean distance of 
        output = (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2)
        return np.sqrt(output);