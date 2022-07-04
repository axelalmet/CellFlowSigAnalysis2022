library(dplyr)
library(CellChat)
library(reticulate)

ad <- import("anndata", convert = FALSE) # Import the module used to read in H5AD files

### Set the ligand-receptor database. Here we will use the "Secreted Signaling" database for cell-cell communication (let's look at ECM-receptor in the future )
CellChatDB <- CellChatDB.mouse # The othe roption is CellChatDB.mouse
showDatabaseCategory(CellChatDB)
CellChatDB.use <- CellChatDB
CellChatDB.subset <- subsetDB(CellChatDB, search = "Secreted Signaling") # Other options include ECM-Receptor and Cell-Cell Contact

min_percentage <- 0.1

### Visium first
visium_ad_object <- ad$read_h5ad("/Users/axelalmet/Documents/CellCellCommunicationModelling/Data/Chen2020/chen20_merged.h5ad")

visium.data.input <- t(py_to_r(visium_ad_object$X))
rownames(visium.data.input) <- rownames(py_to_r(visium_ad_object$var))
colnames(visium.data.input) <- rownames(py_to_r(visium_ad_object$obs))

meta.data <- py_to_r(visium_ad_object$obs)
# Remove NA assignments
meta.data <- meta.data[!is.na(meta.data$AT),]

visium.data.use <- rownames(meta.data)[!is.na(meta.data$AT)] # extract the cell names from disease data
visium.data.input <- visium.data.input[, visium.data.use]

visium.identity <- data.frame(group = meta.data$AT, row.names = row.names(meta.data))

wt.use <- rownames(meta.data)[(meta.data$Group == "WT_03")] # extract the cell names from disease data
ad.use <- rownames(meta.data)[meta.data$Group == "AD_03"] # extract the cell names from disease data

wt.data.input <- visium.data.input[, wt.use]
ad.data.input <- visium.data.input[, ad.use]

wt.meta <- meta.data[meta.data$Group == "WT_03", ]
ad.meta <- meta.data[meta.data$Group == "AD_03", ]

wt.control.identity <- data.frame(group = wt.meta$AT, row.names = row.names(wt.meta))
wt.control.identity$group <- droplevels(wt.control.identity$group)
ad.stim.identity <- data.frame(group = ad.meta$AT, row.names = row.names(ad.meta))
ad.stim.identity$group <- droplevels(ad.stim.identity$group)

# Create the cellchat objects
wt.cc <- createCellChat(object = wt.data.input, do.sparse = T, meta = wt.control.identity, group.by = "group")
levels(wt.cc@idents) # show factor levels of the cell labels

ad.cc <-createCellChat(object = ad.data.input, do.sparse = T, meta = ad.stim.identity, group.by = "group")
levels(ad.cc@idents) # show factor levels of the cell labels

wtGroupSize <- as.numeric(table(wt.cc@idents)) # Get the number of cells in each group
adGroupSize <- as.numeric(table(ad.cc@idents)) # Get the number of cells in each group

# Get the min size to filter out communications
wtMinCells <- floor(0.1 / 100.0 * sum(wtGroupSize))
adMinCells <- floor(0.1 / 100.0 * sum(adGroupSize))

# Set the databases
wt.cc@DB <- CellChatDB.subset 
ad.cc@DB <- CellChatDB.subset 

# We now identify over-expressed ligands/receptors in a cell group and then project gene expression data onto the protein-protein interaction network
### Control
wt.cc <- subsetData(wt.cc) # We subset the expression data of signalling genes to save on computational cost
wt.cc <- identifyOverExpressedGenes(wt.cc) # Identify over-expressed genes (I wonder how much the pre-processing in Seurat has an effect on this)
wt.cc <- identifyOverExpressedInteractions(wt.cc) # Identify the over-expressed ligand-receptor interactions, which are determined by an over-expressed ligand OR receptor
wt.cc <- projectData(wt.cc, PPI.mouse) # Other option includes PPI.human We're told that we may have to comment these out.

wt.cc <- computeCommunProb(wt.cc, raw.use = FALSE, population.size = FALSE)
wt.cc <- filterCommunication(wt.cc, min.cells = wtMinCells) # Filter out the clusters with small numbers
wt.cc <- computeCommunProbPathway(wt.cc) # Calculate the probabilities at the signalling level
wt.cc <- aggregateNet(wt.cc) # Calculates the aggregated cell-cell communication network by counting the links or summing the communication probabilities

### Stimulated
ad.cc <- subsetData(ad.cc) # We subset the expression data of signalling genes to save on computational cost
ad.cc <- identifyOverExpressedGenes(ad.cc) # Identify over-expressed genes (I wonder how much the pre-processing in Seurat has an effect on this)
ad.cc <- identifyOverExpressedInteractions(ad.cc) # Identify the over-expressed ligand-receptor interactions, which are determined by an over-expressed ligand OR receptor
ad.cc <- projectData(ad.cc, PPI.mouse) # Other option includes PPI.human We're told that we may have to comment these out.

ad.cc <- computeCommunProb(ad.cc, raw.use = FALSE, population.size = FALSE)
ad.cc <- filterCommunication(ad.cc, min.cells = adMinCells) # Filter out the clusters with small numbers
ad.cc <- computeCommunProbPathway(ad.cc) # Calculate the probabilities at the signalling level
ad.cc <- aggregateNet(ad.cc) # Calculates the aggregated cell-cell communication network by counting the links or summing the communication probabilities

wt.control.net <- subsetCommunication(wt.cc)
ad.stim.net <- subsetCommunication(ad.cc)

write.csv(wt.control.net, file = "/Users/axelalmet/Documents/scRNASeqAnalysisAndModelling/CellChat/Output/Chen2020/chen20_communications_WT_03.csv")
write.csv(ad.stim.net, file = "/Users/axelalmet/Documents/scRNASeqAnalysisAndModelling/CellChat/Output/Chen2020/chen20_communications_AD_03.csv")
