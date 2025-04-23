------------------------------------------------------------------------------------------------------------------------------------
# ParticleChromo3D -  A Particle Swarm Optimization Algorithm for Chromosome 3D Structure Prediction from Hi-C Data   
------------------------------------------------------------------------------------------------------------------------------------
**OluwadareLab,**
**University of Colorado, Colorado Springs**

----------------------------------------------------------------------
**Developers:** <br />
		 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;David Vadnais<br />
		 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Department of Computer Science <br />
		 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;University of Colorado, Colorado Springs <br />
		 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Email: dvadnais@uccs.edu <br /><br />
		 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Michael Middleton<br />
		 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Department of Computer Science <br />
		 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;University of Colorado, Colorado Springs <br />
		 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Email: mmiddlet@uccs.edu 

**Contact:** <br />
		 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Oluwatosin Oluwadare, PhD <br />
		 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Department of Computer Science <br />
		 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;University of Colorado, Colorado Springs <br />
		 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Email: ooluwada@uccs.edu 
    
--------------------------------------------------------------------	

**0.	Follow-on:**
-----------------------------------------------------------
This code has continued development in the follow on [ParticleChromo3D+](https://github.com/OluwadareLab/ParticleChromo3D_Plus).  As part of that process we released a containerized web gui with restful endpoints that can be accessed [here](http://ParticleChromo3D.online). Additionally, we released a PyPi project which we will describe below.

**ParticleChromo3D from PyPi**

*NOTE: This usage is underconstruction.*<br>
We released a [pypi package.](https://pypi.org/project/ParticleChromo3D/)

Simply run : `pip install ParticleChromo3D`

Then :
```python
from ParticleChromo3D import Ps
import numpy as np

fout = Ps.strip_file("exampleIfs/chr22_matrix.txt")

# this will all get moved into a self contained function soon
theseAlphas = np.array([0.1, 2.0, 0.1]) * 100
theAlphas = ( np.array(range(int(theseAlphas[0]), int(theseAlphas[1]), int(theseAlphas[2]))) / 100)

outputOfSwarm = Ps.Full_List(fout, "this_pdb", theseAlphas)[0]

bestSpearm = outputOfSwarm[1]
bestCost = outputOfSwarm[2]
bestAlpha = theAlphas[outputOfSwarm[4]]
bestPearsonRHO = outputOfSwarm[0]

print(f"Convert factor:: {bestAlpha}")
print(f"SSE at best spearman : {bestCost}")
print(f"Best Spearman correlation Dist vs. Reconstructed Dist  : {bestSpearm}")
print(f"Best Pearson correlation Dist vs. Reconstructed Dist: {bestPearsonRHO}")

Ps.Write_Log(
    "this_run.log", fout, bestAlpha, bestCost, bestSpearm, bestPearsonRHO
)
```

Additional usage documentation can be found in ParticleChromo3D+. [Click here to view](https://github.com/OluwadareLab/ParticleChromo3D_Plus/tree/main/help)

**1.	Content of folders:**
-----------------------------------------------------------	
* src: Python Source Code and utility's used. <br />
* input-and-models: Synthetic and Real Hi-C datasets used. <br />
* Results: Output Structions generated for all the experiments performed.<br />

**2.	Hi-C Data used in this study:**
-----------------------------------------------------------
In our study, we used the synthetic dataset from [Adhikari, et al](https://doi.org/10.1186/s12864-016-3210-4). The contact maps, the original models and their reconstructed models used in this study were downloaded from [here](http://sysbio.rnet.missouri.edu/bdm_download/chromosome3d/unzipped/Input/Synthetic/)

The GM12878 cell Hi-C dataset, GEO Accession number GSE63525, was downloaded from [GSDB](http://sysbio.rnet.missouri.edu/3dgenome/GSDB/details.php?id=GM12878) with GSDB ID: OO7429SF

**3.	Input matrix file format:**
-----------------------------------------------------------

Square Matrix Input format: The square matrix is a space seperated N by N intra-chromosomal contact matrix derived from Hi-C data, where N is the number of equal-sized regions of a chromosome.

**4.	Dependencies Installation:**
-----------------------------------------------------------

Install env with anaconda from spec-file: <br />
conda create --name ParticleChromo3D --file ParticleChromo3d_env.txt

Or install dependencies individually <br />

Dependencies:<br />
&nbsp;&nbsp;&nbsp;&nbsp;Biopython - 1.7.8 <br />
&nbsp;&nbsp;&nbsp;&nbsp;scipy - 1.5.2 <br />
&nbsp;&nbsp;&nbsp;&nbsp;numpy - 1.19.2 <br />

**5.	Usage:**
----------------------------------------------------------- 
usage: **python3 Ps.py infile [-h] [-ss SWARMSIZE] [-itt ITTCOUNT] [-t THRESHOLD] [-rr RANDRANGE] [-lf LOSSFUNCTION]** <br /> 	
                           		
* **positional arguments**: <br />
&nbsp;&nbsp;&nbsp;&nbsp;infile: Input file in the format of a matrix of contacts <br />

* **optional arguments**: <br />	
	&nbsp;&nbsp;&nbsp;&nbsp;-h, --help  show this help message and exit<br />
	&nbsp;&nbsp;&nbsp;&nbsp;-ss SWARMSIZE, --swarmSize SWARMSIZE <br />
		&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Number of particles in system [Default 15] <br />
	&nbsp;&nbsp;&nbsp;&nbsp;-itt ITTCOUNT, --ittCount ITTCOUNT <br />
		&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Maximum number of iterations before stop [Default 30000] <br />
	&nbsp;&nbsp;&nbsp;&nbsp;-t THRESHOLD, --threshold THRESHOLD <br />
		&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Error threshold before stoping [Default 0.000001] <br />	
	&nbsp;&nbsp;&nbsp;&nbsp;-rr RANDRANGE, --randRange RANDRANGE <br />
		&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Range of x,y,z starting coords. Random value bewtween -randRange,randRange [Default 1] <br />
	&nbsp;&nbsp;&nbsp;&nbsp;-o OUTFILE, --outfile OUTFILE <br />
		&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Filename of the output pdb model  [Default ./chr.pdb]<br />
	&nbsp;&nbsp;&nbsp;&nbsp;-lf LOSSFUNCTION, --lossFunction LOSSFUNCTION <br />
		&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;0 = SSE, 1 = MSE, 2 = RMSE, 3 = Huber [Default 2] <br />
 <br />
	
**6.	Output:**
-----------------------------------------------------------
A pdb file and  log file.

**7.	Publication:**
-----------------------------------------------------------
Vadnais, D., Middleton, M. & Oluwadare, O. ParticleChromo3D: a Particle Swarm Optimization algorithm for chromosome 3D structure prediction from Hi-C data. BioData Mining 15, 19 (2022). https://doi.org/10.1186/s13040-022-00305-x
