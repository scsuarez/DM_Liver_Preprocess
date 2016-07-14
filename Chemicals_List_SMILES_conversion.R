# script for obtaining complete chemical list of compounds examined in Drugmatrix Liver
# from chemical list, obtain SMILES for use in Cheminformatics of project - one direction we may go in
# load eset object
## such as eset_liver_pc2_nsFilterNV.RData
## final eset output of DM Liver Preprocess
load("~/Coding Projects/Working Directories/DM Liver Preprocess/Workspaces_And_Objects/eset_liver_pc2_nsFilterNV.RData")

# load bioconductor base
library(Biobase)
# package for query of db matching chemical names and smiles
library(webchem)
library(rpubchem)

# get phenodata DF from expression set object
PhenodataDF <- pData(eset_liver_pc.2_nsFilterNV_rma_qc)
# get list of chemicals from column "chemical"
ChemList <- PhenodataDF[ ,"CHEMICAL"]
# get uniques, as some chemicals are repeated in arrays
uChemList <- unique(ChemList)
# remove controls, i.e. any entry without a chemical identifier
uChemList <- uChemList[uChemList != ""]
# get pubchem chemical ID returns list
ChemCID <- get_cid(uChemList)
# unlist to vector
ChemCID.unlist <- unlist(ChemCID)
# get chemspider ID for chemical list
chemCSID <- get_csid(uChemList, token = "c2d8bc05-d4d6-4322-be7b-c636b26ddc92")
# from pubchem id vector, return smiles currently returning erroneous output - lots of NAs
chemSMILES <- cir_query(ChemCID.unlist, representation = "smiles")


