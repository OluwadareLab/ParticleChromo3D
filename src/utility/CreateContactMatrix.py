import numpy as np
import pandas as pd

inFilePtr = '../../OthersResults/GSE48262_RAW/GSM1173493_cell-1.txt/GSM1173493_cell-1.txt'
outFilePtr = './GSM1173492_Th1_cell1/'
res = 1000000

# Read in file
df = pd.read_table(inFilePtr, dtype={"chrom1": np.character, "coord1": int, "chrom2": np.character, "coord2": int,})
chromList = df['chrom1'].unique()

print(chromList)

for chromNum in chromList:
	# Get specific chromsome contacts only
	contacts = df.loc[df['chrom1'] == chromNum]
	contacts = contacts.loc[contacts['chrom2'] == chromNum]

	#contacts.to_csv(outFilePtr+'chr'+str(chromNum)+'.csv')
	#asd=asdasd

	# Calculate total bead size
	maxVal = contacts[['coord1','coord2']].values.max()
	minVal = contacts[['coord1','coord2']].values.min()

	beadCount = int(np.ceil((maxVal-minVal)/res))
	print(beadCount)
	
	# Create contact matrix filled with zeros
	contMatrix = np.zeros((beadCount,beadCount))

	# Loop through each contact and count it
	counts = {}
	for index, row in contacts.iterrows():
		# Get each bead number from coordinates
		bead1 = int(np.ceil((row['coord1']-minVal)/res))-1
		bead2 = int(np.ceil((row['coord2']-minVal)/res))-1

		# Adds a contact to dictionary
		if (bead1, bead2) in counts:
			counts[(bead1, bead2)] += 1
		elif (bead2, bead1) in counts:
			counts[(bead2, bead1)] += 1
		else:
			counts[(bead1, bead2)] = 1

	# Loops through dictionary and adds contact value to bead1,bead2 and bead2,bead1 index in contact matrix
	for key, val in counts.items():
		contMatrix[key[0]][key[1]] = val
		contMatrix[key[1]][key[0]] = val

	# Save contact matrix
	print(contMatrix.shape)
	np.savetxt(outFilePtr+'chr'+str(chromNum)+'_'+str(int(res/1000))+'1mb.txt', contMatrix, fmt='%1.1f')