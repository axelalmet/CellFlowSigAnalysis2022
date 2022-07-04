# Load the relevant packages
library(dplyr)
library(Seurat)
library(reticulate)

ad <- import("anndata", convert = FALSE) # Import the module used to read in H5AD files

### Set the ligand-receptor database. Here we will use the "Secreted signalling" database for cell-cell communication (let's look at ECM-receptor in the future )
CellChatDB <- CellChatDB.mouse # The othe roption is CellChatDB.mouse
showDatabaseCategory(CellChatDB)
CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling")) # Other options include ECM-Receptor and Cell-Cell Contact

### Marshall et al. (2022) ###
marshall_ad_object <- ad$read_h5ad("/Users/axelalmet/Documents/CellCellCommunicationModelling/Data/Marshall2022/UMOD/marshall22_umod_merged.h5ad")

marshall.data.input.raw <- py_to_r(marshall_ad_object$layers[['counts']])
marshall.data.input.raw <- as.matrix(marshall.data.input.raw)
marshall.data.input.raw <- as(marshall.data.input.raw, "dgCMatrix")
marshall.data.input.raw <- t(marshall.data.input.raw)

marshall.data.input <- t(py_to_r(marshall_ad_object$X))
rownames(marshall.data.input) <- rownames(py_to_r(marshall_ad_object$var))
colnames(marshall.data.input) <- rownames(py_to_r(marshall_ad_object$obs))

# marshall.data.input <- marshall.data.input - min(marshall.data.input) # Need to make sure the data isn't negative (it is in this case)
# access meta data
marshall.meta.data <- py_to_r(marshall_ad_object$obs)
marshall.var.data <- py_to_r(marshall_ad_object$var)
marshall.identity <- data.frame(group = marshall.meta.data$cell_type, row.names = row.names(marshall.meta.data))

WT_samples <- c("WT_01_1a", "WT_01_1b", "WT_01_1c", "WT_01_1d", "WT_01_1e")
KI_samples <- c("KI_01_1a", "KI_01_1b", "KI_01_1c", "KI_01_1d", "KI_01_1e")

for (i in 1:length(WT_samples))
{
  wt_sub_sample <- WT_samples[i]
  ki_sub_sample <- KI_samples[i]
  
  control.use <- rownames(marshall.meta.data)[marshall.meta.data$sub_sample == wt_sub_sample] # extract the cell names from disease data
  stim.use <- rownames(marshall.meta.data)[marshall.meta.data$sub_sample == ki_sub_sample] # extract the cell names from disease data
  
  control.data.input <- marshall.data.input[, control.use]
  stim.data.input <- marshall.data.input[, stim.use]
  
  control.meta <- marshall.meta.data[marshall.meta.data$sub_sample == wt_sub_sample, ]
  stim.meta <- marshall.meta.data[marshall.meta.data$sub_sample == ki_sub_sample, ]
  
  marshall.control.identity <- data.frame(group = control.meta$cell_type, row.names = row.names(control.meta))
  marshall.control.identity$group <- droplevels(marshall.control.identity$group)
  
  marshall.stim.identity <- data.frame(group = stim.meta$cell_type, row.names = row.names(stim.meta))
  marshall.stim.identity$group <- droplevels(marshall.stim.identity$group)
  
  # Create the cellchat objects
  marshall.control.cc <- createCellChat(object = control.data.input, do.sparse = T, meta = marshall.control.identity, group.by = "group")
  levels(marshall.control.cc@idents) # show factor levels of the cell labels
  
  marshall.stim.cc <-createCellChat(object = stim.data.input, do.sparse = T, meta = marshall.stim.identity, group.by = "group")
  levels(marshall.stim.cc@idents) # show factor levels of the cell labels
  
  controlGroupSize <- as.numeric(table(marshall.control.cc@idents)) # Get the number of cells in each group
  stimGroupSize <- as.numeric(table(marshall.stim.cc@idents)) # Get the number of cells in each group
  
  # Set the databases
  marshall.control.cc@DB <- CellChatDB.use 
  marshall.stim.cc@DB <- CellChatDB.use 
  
  # We now identify over-expressed ligands/receptors in a cell group and then project gene expression data onto the protein-protein interaction network
  ### Control
  marshall.control.cc <- subsetData(marshall.control.cc) # We subset the expression data of signalling genes to save on computational cost
  marshall.control.cc <- identifyOverExpressedGenes(marshall.control.cc) # Identify over-expressed genes (I wonder how much the pre-processing in Seurat has an effect on this)
  marshall.control.cc <- identifyOverExpressedInteractions(marshall.control.cc) # Identify the over-expressed ligand-receptor interactions, which are determined by an over-expressed ligand OR receptor
  marshall.control.cc <- projectData(marshall.control.cc, PPI.mouse) # Other option includes PPI.human We're told that we may have to comment these out.
  
  marshall.control.cc <- computeCommunProb(marshall.control.cc, raw.use = FALSE, population.size = FALSE)
  marshall.control.cc <- computeCommunProbPathway(marshall.control.cc) # Calculate the probabilities at the signalling level
  marshall.control.cc <- aggregateNet(marshall.control.cc) # Calculates the aggregated cell-cell communication network by counting the links or summing the communication probabilities
  
  ### Stimulated
  marshall.stim.cc <- subsetData(marshall.stim.cc) # We subset the expression data of signalling genes to save on computational cost
  marshall.stim.cc <- identifyOverExpressedGenes(marshall.stim.cc) # Identify over-expressed genes (I wonder how much the pre-processing in Seurat has an effect on this)
  marshall.stim.cc <- identifyOverExpressedInteractions(marshall.stim.cc) # Identify the over-expressed ligand-receptor interactions, which are determined by an over-expressed ligand OR receptor
  marshall.stim.cc <- projectData(marshall.stim.cc, PPI.mouse) # Other option includes PPI.human We're told that we may have to comment these out.
  
  marshall.stim.cc <- computeCommunProb(marshall.stim.cc, raw.use = FALSE, population.size = FALSE)
  marshall.stim.cc <- computeCommunProbPathway(marshall.stim.cc) # Calculate the probabilities at the signalling level
  marshall.stim.cc <- aggregateNet(marshall.stim.cc) # Calculates the aggregated cell-cell communication network by counting the links or summing the communication probabilities
  
  marshall.control.net <- subsetCommunication(marshall.control.cc)
  marshall.stim.net <- subsetCommunication(marshall.stim.cc)
  
  write.csv(marshall.control.net, file = paste("/Users/axelalmet/Documents/scRNASeqAnalysisAndModelling/CellChat/Output/Marshall2022/UMOD/marshall22_communications_", wt_sub_sample, ".csv", sep=""))
  write.csv(marshall.stim.net, file = paste("/Users/axelalmet/Documents/scRNASeqAnalysisAndModelling/CellChat/Output/Marshall2022/UMOD/marshall22_communications_", ki_sub_sample, ".csv", sep=""))

}