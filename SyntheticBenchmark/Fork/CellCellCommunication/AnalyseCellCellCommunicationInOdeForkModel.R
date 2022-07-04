# Load the relevant packages
library(dplyr)
library(ggplot2)
library(CellChat)
library(reticulate)

ad <- import("anndata", convert = FALSE) # Import the module used to read in H5AD files

### Set the ligand-receptor database. Here we will use the "Secreted signalling" database for cell-cell communication (let's look at ECM-receptor in the future )
CellChatDB <- CellChatDB.human # The othe roption is CellChatDB.human
showDatabaseCategory(CellChatDB)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # Other options include ECM-Receptor and Cell-Cell Contact

fork_ad_object <- ad$read_h5ad("../adata_fork.h5ad")

fork.data.input <- t(py_to_r(fork_ad_object$X))
rownames(fork.data.input) <- as.factor(rownames(py_to_r(fork_ad_object$var)))
colnames(fork.data.input) <- rownames(py_to_r(fork_ad_object$obs))

# lee.data.input <- lee.data.input - min(lee.data.input) # Need to make sure the data isn't negative (it is in this case)
# access meta data
meta.data <- py_to_r(fork_ad_object$obs)
fork.identity <- data.frame(group = meta.data$type, row.names = row.names(meta.data))

obs.use <- rownames(meta.data)[meta.data$Condition == "Observation"] # extract the cell names from disease data
int.use <- rownames(meta.data)[meta.data$Condition == "Intervention"] # extract the cell names from disease data

obs.data.input <- fork.data.input[, obs.use]
int.data.input <- fork.data.input[, int.use]

obs.meta <- meta.data[meta.data$Condition == "Observation", ]
int.meta <- meta.data[meta.data$Condition == "Intervention", ]

obs.identity <- data.frame(group = obs.meta$type, row.names = row.names(obs.meta))
use.identity <- data.frame(group = int.meta$type, row.names = row.names(int.meta))

# Create the cellchat objects
obs.cc <-createCellChat(object = obs.data.input, do.sparse = F, meta = obs.identity, group.by = "group")
levels(obs.cc@idents) # show factor levels of the cell labels

int.cc <-createCellChat(object = int.data.input, do.sparse = F, meta = use.identity, group.by = "group")
levels(int.cc@idents) # show factor levels of the cell labels

obsGroupSize <- as.numeric(table(obs.cc@idents)) # Get the number of cells in each group
intGroupSize <- as.numeric(table(int.cc@idents)) # Get the number of cells in each group

# Set the databases
obs.cc@DB <- CellChatDB.use 
int.cc@DB <- CellChatDB.use 

# We now identify over-expressed ligands/receptors in a cell group and then project gene expression data onto the protein-protein interaction network
### Healthy
obs.cc <- subsetData(obs.cc) # We subset the expression data of signalling genes to save on computational cost
obs.cc <- identifyOverExpressedGenes(obs.cc) # Identify over-expressed genes (I wonder how much the pre-processing in Seurat has an effect on this)
obs.cc <- identifyOverExpressedInteractions(obs.cc) # Identify the over-expressed ligand-receptor interactions, which are determined by an over-expressed ligand OR receptor
# obs.cc <- projectData(obs.cc, PPI.human) # Other option includes PPI.human We're told that we may have to comment these out.

obs.cc <- computeCommunProb(obs.cc, raw.use = TRUE, population.size = FALSE)
obs.cc <- filterCommunication(obs.cc, min.cells = 0.001*sum(obsGroupSize))
obs.cc <- computeCommunProbPathway(obs.cc) # Calculate the probabilities at the signalling level
obs.cc <- aggregateNet(obs.cc) # Calculates the aggregated cell-cell communication network by counting the links or summing the communication probabilities

obs.net <- subsetCommunication(obs.cc)
write.csv(obs.net, file = "../Output/fork_communications_observation.csv")

### Diseased
int.cc <- subsetData(int.cc) # We subset the expression data of signalling genes to save on computational cost
int.cc <- identifyOverExpressedGenes(int.cc) # Identify over-expressed genes (I wonder how much the pre-processing in Seurat has an effect on this)
int.cc <- identifyOverExpressedInteractions(int.cc) # Identify the over-expressed ligand-receptor interactions, which are determined by an over-expressed ligand OR receptor
# int.cc <- projectData(int.cc, PPI.human) # Other option includes PPI.human We're told that we may have to comment these out.

int.cc <- computeCommunProb(int.cc, raw.use = TRUE, population.size = FALSE)
int.cc <- filterCommunication(int.cc, min.cells = 0.001*sum(intGroupSize))
int.cc <- computeCommunProbPathway(int.cc) # Calculate the probabilities at the signalling level
int.cc <- aggregateNet(int.cc) # Calculates the aggregated cell-cell communication network by counting the links or summing the communication probabilities

int.net <- subsetCommunication(int.cc)
write.csv(int.net, file = "../Output/fork_communications_intervention.csv")

