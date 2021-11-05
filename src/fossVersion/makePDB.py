## -- mat2PDB.m --
#
# this function creates a PDB from coordinate data. Represent all inputs as
# a structure field for it to be read. The output format is as given in
# online documentation (as of July 2012 when writing this program)
# http://www.wwpdb.org/documentation/format33/sect9.html#ATOM
#
# Make sure all data input is one-dimensional with the same length. If 
# they are not the same length, the program ignores user input, states 
# an error, and uses defaults. All string inputs should be in cell-format.
# Keep in mind the "element" and "charge" inputs are strings in 
# cell-format, not numbers. 
#
#
# -- required inputs (3) --
#
# input value        meaning
#
# input.X            orthagonal X coordinate data (angstroms)
# input.Y            orthagonal Y coordinate data (angstroms)
# input.Z            orthagonal Z coordinate data (angstroms)
#
# -- optional inputs (12): generates defaults when not user-specified --
#
# input value        meaning                           default value
#
# input.outfile      output file name                 "mat2PDB.pdb"
# input.recordName   output record name of atoms      "ATOM"
# input.atomNum      atom serial number                sequential number
# input.atomName     name of atoms                    "OW" (water oxygen)
# input.altLoc       alt. location indicator          " "
# input.resName      name of residue                  "SOL" (water)
# 
# input.chainID      protein chain identifier         "A"
# input.resNum       residue sequence number           sequential number
# input.occupancy    occupancy factor                 "1.00"
# input.betaFactor   beta factor, temperature         "0.00"
# input.element      element symbol                   "O" (oxygen)
# input.charge       atomic charge                    " "
#
#
# -- example uses --
#
# # translates both X and Y coordinates of 3IJU.pdb by 5 angstroms
# PDBdata = pdb2mat('3IJU.pdb');
# PDBdata.X = PDBdata.X + 5;
# PDBdata.Y = PDBdata.Y + 5;
# mat2pdb(PDBdata)
# 
# # make a PDB with 30 atoms in random places within a 10 angstrom box
# data.X = rand(1,20)*10;
# data.Y = rand(1,20)*10;
# data.Z = rand(1,20)*10;
# mat2pdb(data)
#
# 
import numpy as np

class makePDB:
    def __init__(self):
        ''' Constructor for this class. '''
        
    def fillWithSpace(self, val,sizeGoal,left=False):
        if left == False:
            while len(val) < sizeGoal:
                val = " " + val
        else:
            while len(val) < sizeGoal:
                val =  val + " " 
        return val

    def xyzOutFormat(self, P):
        nP = format(P, '8.3f');
        nP = self.fillWithSpace(nP,8)
        return nP

    def occupancyOutFormat(self, P):
        nP = format(P, '6.2f');
        return nP

    def mat2pdb(self, inputO):
    ## review XYZ coordinate data 
    # coordinate data is required! Checking XYZ input
        if len(inputO.xyzData[0,:]) != 3:
            fprintf('we need xyz coordinate data to make a PDB!!\n\texiting...\n');
            return;

        X = inputO.xyzData[:,0];
        Y = inputO.xyzData[:,1];
        Z = inputO.xyzData[:,2];


        # in case optional data data not given, fill in blanks
        if ~hasattr(inputO, 'atomNum'):
            inputO.atomNum = range(len(X));

        if ~hasattr(inputO, 'atomName'):
            inputO.atomName = []
            for i in range(len(X)):
                inputO.atomName.append("CA ")

        if ~hasattr(inputO, 'altLoc'):
            inputO.altLoc = []
            for i in range(len(X)):
                inputO.altLoc.append(" ")

        if ~hasattr(inputO, 'resName'):
            inputO.resName = []
            for i in range(len(X)):
                inputO.resName.append("MET ")

        if ~hasattr(inputO, 'chainID'):
            inputO.chainID = []
            for i in range(len(X)):
                inputO.chainID.append('B')

        if ~hasattr(inputO, 'resNum'):
            inputO.resNum = range(len(X));

        if ~hasattr(inputO, 'occupancy'):
            inputO.occupancy = np.ones(len(X));


        if ~hasattr(inputO, 'betaFactor'):
            inputO.betaFactor = np.zeros(len(X));

        if ~hasattr(inputO, 'element'):
            inputO.element = []
            for i in range(len(X)):
                inputO.element.append('O')

        if ~hasattr(inputO, 'charge'):
            inputO.charge = []
            for i in range(len(X)):
                inputO.charge.append(" ")

        outfile    = inputO.outputFile;
        recordName = inputO.recordName;
        atomNum    = inputO.atomNum;
        atomName   = inputO.atomName;
        altLoc     = inputO.altLoc;
        resName    = inputO.resName;
        chainID    = inputO.chainID;
        resNum     = inputO.resNum;
        occupancy  = inputO.occupancy;
        betaFactor = inputO.betaFactor;## temperature

        element    = inputO.element;
        charge     = inputO.charge;

        ## remove faulty inputs
            ### TODO

        ## create PDB

        # open file
        print('outputting PDB in file : ', outfile);
        FILE = open(outfile, 'w');
        aNumStrCorrected = []
        resNumStrCorrected = []

        # output data
        for n in range(len(atomNum)):
            # standard PDB output line


            tmpStr = str(resNum[n])
            while len(tmpStr) < 3:
                tmpStr = " " + tmpStr
            resNumStrCorrected.append(tmpStr)


            inString = (recordName[n] + "  " + self.fillWithSpace(str(atomNum[n]),5) +" "
                + self.fillWithSpace(atomName[n],4,True) + 
                altLoc[n] + resName[n] + chainID[n] +
                self.fillWithSpace(str(resNum[n]),4) + " "#Code for insertions of residues 27
                +"   " +self.xyzOutFormat(X[n]) + 
                self.xyzOutFormat(Y[n])+ self.xyzOutFormat(Z[n])
                + self.fillWithSpace(self.occupancyOutFormat(occupancy[n]),6) 
                + self.fillWithSpace(self.occupancyOutFormat(betaFactor[n]),6) +"       "#7 spaces
                + "    " #Segment identifier
                + self.fillWithSpace(element[n],3, True ) 
                + self.fillWithSpace(charge[n],2) +'\n'
            )



            FILE.write(inString);

            # output progress in terminal
            #if ~(n % 400):
            #    print('Progress : ', 100*n / len(atomNum));

        for n in range(len(atomNum)-1):
            inString2 = 'CONECT' + self.fillWithSpace(str(atomNum[n]),5) + self.fillWithSpace(str(atomNum[n+1]),5)+ '\n'
            FILE.write(inString2);

        FILE.write('END\n');

        # close file
        print('done! closing file...\n');

        FILE.close();

