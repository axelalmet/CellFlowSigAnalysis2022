# Load the relevant packages
library(dplyr)
library(CellChat)
library(reticulate)

ad <- import("anndata", convert = FALSE) # Import the module used to read in H5AD files

### Set the ligand-receptor database. Here we will use the "Secreted signalling" database for cell-cell communication (let's look at ECM-receptor in the future )
CellChatDB <- CellChatDB.human # The othe roption is CellChatDB.human
showDatabaseCategory(CellChatDB)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # Other options include ECM-Receptor and Cell-Cell Contact

### Kang et al. (2018) ###
kang_ad_object <- ad$read_h5ad("/Users/axelalmet/Documents/CellCellCommunicationModelling/Data/Kang2018/kang18_tutorial.h5ad")

kang.data.input <- t(py_to_r(kang_ad_object$X))
rownames(kang.data.input) <- rownames(py_to_r(kang_ad_object$var))
colnames(kang.data.input) <- rownames(py_to_r(kang_ad_object$obs))

# kang.data.input <- kang.data.input - min(kang.data.input) # Need to make sure the data isn't negative (it is in this case)
# access meta data
meta.data <- py_to_r(kang_ad_object$obs)
kang.identity <- data.frame(group = meta.data$cell_type, row.names = row.names(meta.data))

control.use <- rownames(meta.data)[meta.data$condition == "control"] # extract the cell names from disease data
stim.use <- rownames(meta.data)[meta.data$condition == "stimulated"] # extract the cell names from disease data

control.data.input <- kang.data.input[, control.use]
stim.data.input <- kang.data.input[, stim.use]

control.meta <- meta.data[meta.data$condition == "control", ]
stim.meta <- meta.data[meta.data$condition == "stimulated", ]

kang.control.identity <- data.frame(group = control.meta$cell_type, row.names = row.names(control.meta))
kang.stim.identity <- data.frame(group = stim.meta$cell_type, row.names = row.names(stim.meta))

# Create the cellchat objects
kang.control.cc <-createCellChat(object = control.data.input, do.sparse = T, meta = kang.control.identity, group.by = "group")
levels(kang.control.cc@idents) # show factor levels of the cell labels

kang.stim.cc <-createCellChat(object = stim.data.input, do.sparse = T, meta = kang.stim.identity, group.by = "group")
levels(kang.stim.cc@idents) # show factor levels of the cell labels

controlGroupSize <- as.numeric(table(kang.control.cc@idents)) # Get the number of cells in each group
stimGroupSize <- as.numeric(table(kang.stim.cc@idents)) # Get the number of cells in each group

# Set the databases
kang.control.cc@DB <- CellChatDB.use 
kang.stim.cc@DB <- CellChatDB.use 

# We now identify over-expressed ligands/receptors in a cell group and then project gene expression data onto the protein-protein interaction network
### Control
kang.control.cc <- subsetData(kang.control.cc) # We subset the expression data of signalling genes to save on computational cost
kang.control.cc <- identifyOverExpressedGenes(kang.control.cc) # Identify over-expressed genes (I wonder how much the pre-processing in Seurat has an effect on this)
kang.control.cc <- identifyOverExpressedInteractions(kang.control.cc) # Identify the over-expressed ligand-receptor interactions, which are determined by an over-expressed ligand OR receptor
kang.control.cc <- projectData(kang.control.cc, PPI.human) # Other option includes PPI.human We're told that we may have to comment these out.

kang.control.cc <- computeCommunProb(kang.control.cc, raw.use = FALSE, population.size = FALSE)
kang.control.cc <- computeCommunProbPathway(kang.control.cc) # Calculate the probabilities at the signalling level
kang.control.cc <- aggregateNet(kang.control.cc) # Calculates the aggregated cell-cell communication network by counting the links or summing the communication probabilities

### Stimulated
kang.stim.cc <- subsetData(kang.stim.cc) # We subset the expression data of signalling genes to save on computational cost
kang.stim.cc <- identifyOverExpressedGenes(kang.stim.cc) # Identify over-expressed genes (I wonder how much the pre-processing in Seurat has an effect on this)
kang.stim.cc <- identifyOverExpressedInteractions(kang.stim.cc) # Identify the over-expressed ligand-receptor interactions, which are determined by an over-expressed ligand OR receptor
kang.stim.cc <- projectData(kang.stim.cc, PPI.human) # Other option includes PPI.human We're told that we may have to comment these out.

kang.stim.cc <- computeCommunProb(kang.stim.cc, raw.use = FALSE, population.size = FALSE)
kang.stim.cc <- computeCommunProbPathway(kang.stim.cc) # Calculate the probabilities at the signalling level
kang.stim.cc <- aggregateNet(kang.stim.cc) # Calculates the aggregated cell-cell communication network by counting the links or summing the communication probabilities

kang.control.net <- subsetCommunication(kang.control.cc)
kang.stim.net <- subsetCommunication(kang.stim.cc)

write.csv(kang.control.net, file = "/Users/axelalmet/Documents/scRNASeqAnalysisAndModelling/CellChat/Output/Kang2018/kang18_communications_control.csv")
write.csv(kang.stim.net, file = "/Users/axelalmet/Documents/scRNASeqAnalysisAndModelling/CellChat/Output/Kang2018/kang18_communications_stimulated.csv")
