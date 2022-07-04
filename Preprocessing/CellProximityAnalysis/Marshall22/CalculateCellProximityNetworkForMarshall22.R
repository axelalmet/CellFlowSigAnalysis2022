library(Giotto)
library(reticulate)
library(Matrix)

ad <- import("anndata", convert = FALSE) # Import the module used to read in H5AD files

### Marshall et al. (2022) ###
marshall_data_directory <- "/Users/axelalmet/Documents/CellCellCommunicationModelling/Data/Marshall2022/UMOD/"
marshall_ad_object <- ad$read_h5ad("/Users/axelalmet/Documents/CellCellCommunicationModelling/Data/Marshall2022/UMOD/marshall22_umod_merged.h5ad")

marshall.data.input.raw <- py_to_r(marshall_ad_object$layers[['counts']])
marshall.data.input.raw <- as.matrix(marshall.data.input.raw)
marshall.data.input.raw <- as(marshall.data.input.raw, "dgCMatrix")
marshall.data.input.raw <- t(marshall.data.input.raw)

marshall.data.input <- py_to_r(marshall_ad_object$X)
marshall.data.input <- as.matrix(marshall.data.input)
marshall.data.input <- as(marshall.data.input, "dgCMatrix")
marshall.data.input <- t(marshall.data.input)

marshall.meta.data <- py_to_r(marshall_ad_object$obs)
marshall.var.data <- py_to_r(marshall_ad_object$var)
# Remove NA assignments
marshall.meta.data <- marshall.meta.data[!is.na(marshall.meta.data$cell_type),]

rownames(marshall.data.input) <- rownames(var.data)
colnames(marshall.data.input) <- rownames(marshall.meta.data)
rownames(marshall.data.input.raw) <- rownames(var.data)
colnames(marshall.data.input.raw) <- rownames(marshall.meta.data)

marshall.data.use <- rownames(marshall.meta.data)[!is.na(marshall.meta.data$cell_type)] # extract the cell names from disease data
marshall.data.input <- marshall.data.input[, marshall.data.use]
marshall.data.input.raw <- marshall.data.input.raw[, marshall.data.use]

WT_samples <- c("WT_01_1a", "WT_01_1b", "WT_01_1c", "WT_01_1d", "WT_01_1e")
KI_samples <- c("KI_01_1a", "KI_01_1b", "KI_01_1c", "KI_01_1d", "KI_01_1e")

my_instructions = createGiottoInstructions(python_path = '/usr/local/Cellar/python@3.8/3.8.13/bin/python3')

for (i in 1:length(WT_samples))
{
  wt_sub_sample <- WT_samples[i]
  ki_sub_sample <- KI_samples[i]
  
  wt.use <- rownames(marshall.meta.data)[(marshall.meta.data$sub_sample == wt_sub_sample)] # extract the cell names from disease data
  ki.use <- rownames(marshall.meta.data)[marshall.meta.data$sub_sample == ki_sub_sample] # extract the cell names from disease data
  
  wt.data.input <- marshall.data.input[, wt.use]
  ki.data.input <- marshall.data.input[, ki.use]
  
  wt.meta <- marshall.meta.data[marshall.meta.data$sub_sample == wt_sub_sample, ]
  ki.meta <- marshall.meta.data[marshall.meta.data$sub_sample == ki_sub_sample, ]
  
  wt.spatial.positions <- cbind(wt.meta$xcoord, wt.meta$ycoord)
  rownames(wt.spatial.positions) <- wt.use
  ki.spatial.positions <- cbind(ki.meta$xcoord, ki.meta$ycoord)
  rownames(ki.spatial.positions) <- ki.use
  
  wt_giotto_object = createGiottoObject(raw_exprs = wt.data.input,
                                        spatial_locs = wt.spatial.positions,
                                        cell_metadata = wt.meta,
                                        gene_metadata = marshall.var.data,
                                        instructions = my_instructions)
  
  ki_giotto_object = createGiottoObject(raw_exprs = ki.data.input,
                                        spatial_locs = ki.spatial.positions,
                                        cell_metadata = ki.meta,
                                        gene_metadata = marshall.var.data,
                                        instructions = my_instructions)
  
  # Plot the cell proximity network for the WT data
  plotStatDelaunayNetwork(gobject = wt_giotto_object, maximum_distance = 400, save_plot = F)
  wt_giotto_object = createSpatialNetwork(gobject = wt_giotto_object, minimum_k = 2, maximum_distance_delaunay = 400)
  
  wt_cell_proximities = cellProximityEnrichment(gobject = wt_giotto_object,
                                                cluster_column = 'cell_type',
                                                adjust_method = 'fdr',
                                                number_of_simulations = 2000)
  
  
  # Plot the cell proximity network for the AD data
  plotStatDelaunayNetwork(gobject = ki_giotto_object, maximum_distance = 400, save_plot = F)
  ki_giotto_object = createSpatialNetwork(gobject = ki_giotto_object, minimum_k = 2, maximum_distance_delaunay = 400)
  
  ki_cell_proximities = cellProximityEnrichment(gobject = ki_giotto_object,
                                                cluster_column = 'cell_type',
                                                adjust_method = 'fdr',
                                                number_of_simulations = 2000)

  
  # Save the files
  write.csv(wt_cell_proximities$enrichm_res, paste(marshall_data_directory, "marshall22_cellproximities_", wt_sub_sample, ".csv", sep=""))
  write.csv(ki_cell_proximities$enrichm_res, paste(marshall_data_directory, "marshall22_cellproximities_", ki_sub_sample, ".csv", sep=""))
  
}
