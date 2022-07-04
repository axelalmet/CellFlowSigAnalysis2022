# Load the relevant packages
library(dplyr)
library(CellChat)
library(reticulate)

ad <- import("anndata", convert = FALSE) # Import the module used to read in H5AD files

### Set the ligand-receptor database. Here we will use the "Secreted signalling" database for cell-cell communication (let's look at ECM-receptor in the future )
CellChatDB <- CellChatDB.human # The othe roption is CellChatDB.human
showDatabaseCategory(CellChatDB)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # Other options include ECM-Receptor and Cell-Cell Contact

### mcfaline et al. (2018) ###
mcfaline_ad_object <- ad$read_h5ad("/Users/axelalmet/Documents/CellCellCommunicationModelling/Data/McFaline-Figueroa2019/mcfaline19_merged.h5ad")

mcfaline.data.input <- t(py_to_r(mcfaline_ad_object$X))
rownames(mcfaline.data.input) <- rownames(py_to_r(mcfaline_ad_object$var))
colnames(mcfaline.data.input) <- rownames(py_to_r(mcfaline_ad_object$obs))

# mcfaline.data.input <- mcfaline.data.input - min(mcfaline.data.input) # Need to make sure the data isn't negative (it is in this case)
# access meta data
meta.data <- py_to_r(mcfaline_ad_object$obs)
mcfaline.identity <- data.frame(group = meta.data$leiden, row.names = row.names(meta.data))

inner.use <- rownames(meta.data)[meta.data$spatial_id == "inner"] # extract the cell names from disease data
outer.use <- rownames(meta.data)[meta.data$spatial_id == "outer"] # extract the cell names from disease data

inner.data.input <- mcfaline.data.input[, inner.use]
outer.data.input <- mcfaline.data.input[, outer.use]

inner.meta <- meta.data[meta.data$spatial_id == "inner", ]
outer.meta <- meta.data[meta.data$spatial_id == "outer", ]

mcfaline.inner.identity <- data.frame(group = inner.meta$leiden, row.names = row.names(inner.meta))
mcfaline.outer.identity <- data.frame(group = outer.meta$leiden, row.names = row.names(outer.meta))

# Create the cellchat objects
mcfaline.inner.cc <-createCellChat(object = inner.data.input, do.sparse = T, meta = mcfaline.inner.identity, group.by = "group")
levels(mcfaline.inner.cc@idents) # show factor levels of the cell labels

mcfaline.outer.cc <-createCellChat(object = outer.data.input, do.sparse = T, meta = mcfaline.outer.identity, group.by = "group")
levels(mcfaline.outer.cc@idents) # show factor levels of the cell labels

innerGroupSize <- as.numeric(table(mcfaline.inner.cc@idents)) # Get the number of cells in each group
stimGroupSize <- as.numeric(table(mcfaline.outer.cc@idents)) # Get the number of cells in each group

# Set the databases
mcfaline.inner.cc@DB <- CellChatDB.use 
mcfaline.outer.cc@DB <- CellChatDB.use 

# We now identify over-expressed ligands/receptors in a cell group and then project gene expression data onto the protein-protein interaction network
### inner
mcfaline.inner.cc <- subsetData(mcfaline.inner.cc) # We subset the expression data of signalling genes to save on computational cost
mcfaline.inner.cc <- identifyOverExpressedGenes(mcfaline.inner.cc) # Identify over-expressed genes (I wonder how much the pre-processing in Seurat has an effect on this)
mcfaline.inner.cc <- identifyOverExpressedInteractions(mcfaline.inner.cc) # Identify the over-expressed ligand-receptor interactions, which are determined by an over-expressed ligand OR receptor
mcfaline.inner.cc <- projectData(mcfaline.inner.cc, PPI.human) # Other option includes PPI.human We're told that we may have to comment these out.

mcfaline.inner.cc <- computeCommunProb(mcfaline.inner.cc, raw.use = FALSE, population.size = FALSE)
mcfaline.inner.cc <- computeCommunProbPathway(mcfaline.inner.cc) # Calculate the probabilities at the signalling level
mcfaline.inner.cc <- aggregateNet(mcfaline.inner.cc) # Calculates the aggregated cell-cell communication network by counting the links or summing the communication probabilities

### outer
mcfaline.outer.cc <- subsetData(mcfaline.outer.cc) # We subset the expression data of signalling genes to save on computational cost
mcfaline.outer.cc <- identifyOverExpressedGenes(mcfaline.outer.cc) # Identify over-expressed genes (I wonder how much the pre-processing in Seurat has an effect on this)
mcfaline.outer.cc <- identifyOverExpressedInteractions(mcfaline.outer.cc) # Identify the over-expressed ligand-receptor interactions, which are determined by an over-expressed ligand OR receptor
mcfaline.outer.cc <- projectData(mcfaline.outer.cc, PPI.human) # Other option includes PPI.human We're told that we may have to comment these out.

mcfaline.outer.cc <- computeCommunProb(mcfaline.outer.cc, raw.use = FALSE, population.size = FALSE)
mcfaline.outer.cc <- computeCommunProbPathway(mcfaline.outer.cc) # Calculate the probabilities at the signalling level
mcfaline.outer.cc <- aggregateNet(mcfaline.outer.cc) # Calculates the aggregated cell-cell communication network by counting the links or summing the communication probabilities

mcfaline.inner.net <- subsetCommunication(mcfaline.inner.cc)
mcfaline.outer.net <- subsetCommunication(mcfaline.outer.cc)

write.csv(mcfaline.inner.net, file = "/Users/axelalmet/Documents/scRNASeqAnalysisAndModelling/CellChat/Output/McFaline-Figueroa2019/mcfaline19_communications_inner.csv")
write.csv(mcfaline.outer.net, file = "/Users/axelalmet/Documents/scRNASeqAnalysisAndModelling/CellChat/Output/McFaline-Figueroa2019/mcfaline19_communications_outer.csv")
