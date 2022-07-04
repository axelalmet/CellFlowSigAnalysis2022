library(Giotto)
library(reticulate)
library(Matrix)

ad <- import("anndata", convert = FALSE) # Import the module used to read in H5AD files

### Visium first
chen_data_directory <- "/Users/axelalmet/Documents/CellCellCommunicationModelling/Data/Chen2020/"
visium_ad_object <- ad$read_h5ad("/Users/axelalmet/Documents/CellCellCommunicationModelling/Data/Chen2020/chen20_merged.h5ad")

visium.data.input.raw <- py_to_r(visium_ad_object$layers[['counts']])
visium.data.input.raw <- as.matrix(visium.data.input.raw)
visium.data.input.raw <- as(visium.data.input.raw, "dgCMatrix")
visium.data.input.raw <- t(visium.data.input.raw)

visium.data.input <- py_to_r(visium_ad_object$X)
visium.data.input <- as.matrix(visium.data.input)
visium.data.input <- as(visium.data.input, "dgCMatrix")
visium.data.input <- t(visium.data.input)

meta.data <- py_to_r(visium_ad_object$obs)
var.data <- py_to_r(visium_ad_object$var)
# Remove NA assignments
meta.data <- meta.data[!is.na(meta.data$AT),]

rownames(visium.data.input) <- rownames(var.data)
colnames(visium.data.input) <- rownames(meta.data)
rownames(visium.data.input.raw) <- rownames(var.data)
colnames(visium.data.input.raw) <- rownames(meta.data)

visium.data.use <- rownames(meta.data)[!is.na(meta.data$AT)] # extract the cell names from disease data
visium.data.input <- visium.data.input[, visium.data.use]
visium.data.input.raw <- visium.data.input.raw[, visium.data.use]

wt.use <- rownames(meta.data)[(meta.data$Group == "WT_18")] # extract the cell names from disease data
ad.use <- rownames(meta.data)[meta.data$Group == "AD_18"] # extract the cell names from disease data

wt.data.input.raw <- visium.data.input.raw[, wt.use]
ad.data.input.raw <- visium.data.input.raw[, ad.use]

wt.data.input <- visium.data.input[, wt.use]
ad.data.input <- visium.data.input[, ad.use]

wt.meta <- meta.data[meta.data$Group == "WT_18", ]
ad.meta <- meta.data[meta.data$Group == "AD_18", ]

wt.spatial.positions <- cbind(wt.meta$coord_X, wt.meta$coord_Y)
rownames(wt.spatial.positions) <- wt.use
ad.spatial.positions <- cbind(ad.meta$coord_X, ad.meta$coord_Y)
rownames(ad.spatial.positions) <- ad.use

my_instructions = createGiottoInstructions(python_path = '/usr/local/Cellar/python@3.8/3.8.13/bin/python3')

wt_giotto_object = createGiottoObject(raw_exprs = wt.data.input.raw,
                                      spatial_locs = wt.spatial.positions,
                                      cell_metadata = wt.meta,
                                      gene_metadata = var.data,
                                      instructions = my_instructions)

ad_giotto_object = createGiottoObject(raw_exprs = ad.data.input,
                                      spatial_locs = ad.spatial.positions,
                                      cell_metadata = ad.meta,
                                      gene_metadata = var.data,
                                      instructions = my_instructions)

# Plot the cell proximity network for the WT data
plotStatDelaunayNetwork(gobject = wt_giotto_object, maximum_distance = 400, save_plot = F)
wt_giotto_object = createSpatialNetwork(gobject = wt_giotto_object, minimum_k = 2, maximum_distance_delaunay = 400)

wt_cell_proximities = cellProximityEnrichment(gobject = wt_giotto_object,
                                              cluster_column = 'AT',
                                              adjust_method = 'fdr',
                                              number_of_simulations = 2000)

wt_proximity_network <- cellProximityNetwork(gobject = wt_giotto_object, CPscore = wt_cell_proximities,
                                             remove_self_edges = T,
                                             self_loop_strength = 0.2,
                                             only_show_enrichment_edges = F,
                                             rescale_edge_weights = T)

wt_cell_proximity_heatmap <- cellProximityHeatmap(gobject = wt_giotto_object, CPscore = wt_cell_proximities, order_cell_types = F, scale = T,
                                                  color_breaks = c(-2, 0, 1), color_names = c('blue', 'white', 'red'))

# Plot the cell proximity network for the AD data
plotStatDelaunayNetwork(gobject = ad_giotto_object, maximum_distance = 400, save_plot = F)
ad_giotto_object = createSpatialNetwork(gobject = ad_giotto_object, minimum_k = 2, maximum_distance_delaunay = 400)

ad_cell_proximities = cellProximityEnrichment(gobject = ad_giotto_object,
                                              cluster_column = 'AT',
                                              adjust_method = 'fdr',
                                              number_of_simulations = 2000)

ad_proximity_network <- cellProximityNetwork(gobject = ad_giotto_object, CPscore = ad_cell_proximities,
                                             remove_self_edges = T,
                                             self_loop_strength = 0.2,
                                             only_show_enrichment_edges = F,
                                             rescale_edge_weights = T)


ad_cell_proximity_heatmap <- cellProximityHeatmap(gobject = ad_giotto_object, CPscore = ad_cell_proximities, order_cell_types = F, scale = T,
                                                  color_breaks = c(-1, 0, 1), color_names = c('blue', 'white', 'red'))

# Save the files
write.csv(wt_cell_proximities$enrichm_res, paste(chen_data_directory, "chen20_cellproximities_WT_18.csv", sep=""))
write.csv(ad_cell_proximities$enrichm_res, paste(chen_data_directory, "chen20_cellproximities_AD_18.csv", sep=""))

