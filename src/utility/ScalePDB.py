import Helper


inFilePtr = "./output/scale/chrX_1mb_matrix5000.pdb"
compFilePtr = "./output/scale/Th1_Cell1_chrX_model.pdb"

newXYZCommp, newXYZIn = Helper.Proc_PDB(inFilePtr, compFilePtr)

Helper.Write_Output(compFilePtr+'scaled', newXYZCommp)
Helper.Write_Output(inFilePtr+'scaled', newXYZIn)