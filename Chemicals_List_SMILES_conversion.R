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

# get phenodata DF from expression set object using Bioconductor Base package
phenodataDF <- pData(eset_liver_pc.2_nsFilterNV_rma_qc)
# get list of chemicals from column "chemical"
chemList <- phenodataDF[ ,"CHEMICAL"]
# get uniques, as some chemicals are repeated in arrays
uChemList <- unique(chemList)
# remove controls, i.e. any entry without a chemical identifier
uChemList <- uChemList[uChemList != ""]
# get pubchem chemical ID returns list
# chemCID <- get_cid(uChemList)
chemCSID <- cts_convert(uChemList, 'Chemical Name', 'ChemSpider')
chemPC_CID <- cts_convert(uChemList, 'Chemical Name', 'PubChem CID')
chemPC_CID.null <- 
# unlist to vector
chemCSIDdf <- as.data.frame(unlist(chemCSID))
# get chemspider ID for chemical list, token belongs to SCS for chemspider.com
chemCSID <- get_csid(uChemList, token = "c2d8bc05-d4d6-4322-be7b-c636b26ddc92")
# from pubchem id vector, return smiles currently returning erroneous output - lots of NAs
chemSMILES <- cir_query(ChemCID.unlist, representation = "smiles")


