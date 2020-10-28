require(lidR)
require(TreeLS)
require(FNN)
require(RANN)
require(sf)
require(cppRouting)


##------------------------------------------------------------------------------------------------------
##--------------------- GRAPH SEGMENTATION FUNCTION ----------------------------------------------------
##------------------------------------------------------------------------------------------------------
segment.graph <- function(las = NULL, tree.locations = NULL, k = NULL, distance.threshold = 0.38,
                          use.metabolic.scale = FALSE, subsample.graph = 0.2,
                          return.dense = FALSE){
  ## pcd
  segmented <- las

  ## if subsampling, do that now
  if(!is.null(subsample.graph)){
    segmented <- tlsSample(segmented, smp.voxelize(spacing = subsample.graph))
  }
  ## get rid of ground points
  segmented <- filter_poi(segmented, Classification != 2, Z > (1/3))

  ## how many tree objects
  ntree <- nrow(tree.locations)

  ## Give tree locations a uniform Z coordinate
  tree.locations$Z <- 1.0
  tree.locations$Radius[is.na(tree.locations$Radius)] = 0.01
  ## Make sure tree location is close to bole points (if in center, distance threshold could cut it off)
  ## Eventually, this should be a value tied to the actual DBH of the tree
  for(d in 1:nrow(tree.locations)){
    if(tree.locations$Radius[d] > distance.threshold){
      if(tree.locations$Y[d] > mean(tree.locations$Y)){
        tree.locations$Y <- tree.locations$Y - (tree.locations$Radius)
      } else {
        tree.locations$Y <- tree.locations$Y + (tree.locations$Radius)
      }
    }}

  ## ---- Determining metabolic scaling factors for individual trees ----
  ## This is the metabolic scaling factor which will get multiplied by the weights so that smaller trees don't "suck up" points from nearby larger trees
  # met_scale <- 1 / ( tree.locations$Radius^(2/3) )## it seems like the equation in the original manuscript doesn't work that well for me... this can be modified
  met_scale <- (1 / ( tree.locations$Radius))^(3/2) ## it seems like the equation in the original manuscript doesn't work that well for me... this can be modified

  ##---------------------- Identify neighbors  -------------------------------------
  ## create the network by connecting lidar points within X meters
  ## and add in the points for the tree centers to make connections
  point_density <- density(segmented)

  # point_dist_threshold <- round((1/sqrt(point_density))*k, 1) ## in meters, the distance within which points will be connected into a network
  # point_dist_threshold <- (1/3)
  # point_dist_threshold <- 0.42

  ## k is arbitrary... not sure if this should be based on a calculation, e.g. ...
  # k = ((sqrt(point_density))/(1/3))

  ## in meters, the distance within which points will be connected into a network
  k_tree <- get.knn(rbind(segmented@data[ ,c("X","Y","Z")], tree.locations[ ,c("X","Y","Z")]), k = k,  algorithm="cover_tree")

  ##---------------------- Build the Graph OBJ  -------------------------------------
  ## matrix of from, to, weight (dist)
  weight_matrix = cbind(rep(1:dim(k_tree$nn.index)[1],each=k),
                        as.vector(t(k_tree$nn.index)),
                        as.vector(t(k_tree$nn.dist)))
  # which(weight_matrix[,3] > point_dist_threshold)
  # hist(weight_matrix[which(weight_matrix[,3] > point_dist_threshold),3])

  ## get rid of rows in the matrix that have a distance greater than the point_dist_threshold
  weight_matrix <- weight_matrix[-(which(weight_matrix[,3] > distance.threshold)), ]

  ## actually make the graph using this matrix
  graph_obj <- cppRouting::makegraph(weight_matrix,
                                     directed = FALSE) ## a new package

  ## ---- IDs for the origin tree locations ----
  treeloc_ids <- (nrow(segmented@data)+1):(nrow(k_tree$nn.index))
  treeloc_ids <- treeloc_ids[treeloc_ids %in% c(weight_matrix[,1],weight_matrix[,2])]

  ## ---- IDs for the origin tree locations ----
  remaining_ids <- unique(c(weight_matrix[,1],weight_matrix[,2])[!(c(weight_matrix[,1],weight_matrix[,2]) %in% treeloc_ids)])

  ## ---- IDs not used in graph ----
  unused_ids <- (1:(nrow(segmented@data)))[!(1:(nrow(segmented@data)) %in% c(weight_matrix[,1],weight_matrix[,2]))]

  ##-------------- Calculate the shortest paths from all tree origins  -----------------------------
  shortpath <- cppRouting::get_distance_matrix(graph_obj, from = treeloc_ids, to = remaining_ids, allcores = FALSE)
  # shortpath[ ,1:100] ## testing

  # shortpath <- cppRouting::get_multi_paths(graph_obj, from = treeloc_ids, to = remaining_ids, long = TRUE)

  ##-------------- Determine which tree ID is linked to the shortest path  -----------------------------
  ## Use the metabolic scaling factors to modify the distances for each treeID
  if(use.metabolic.scale == TRUE){
    shortpath <- shortpath * met_scale
  }

  ## add this so you don't drop columns in the colMins function
  shortpath <- rbind(shortpath, 999)

  ## Classify the vector of treeIDs (we need all of the point ids back together)
  treeID <- c( Rfast::colMins(shortpath), rep(0, length(unused_ids)) )
  ## give the values the proper point ID name
  names(treeID) <- c(colnames(shortpath), unused_ids)
  ## now get the treeID values back in order to attach to the las object
  treeID <- treeID[order(as.numeric(names(treeID)))]
  ## change the "extra" treeID to 0
  treeID[treeID > ntree] <- 0

  ## -------------------------- add the treeID data to the las object -----------------------------------
  segmented <- add_lasattribute(segmented, (treeID), name = "treeID", desc = "tree id")
  # segmented <- add_lasattribute(segmented, (as.numeric(names(treeID))), name = "name", desc = "tree name")

  ## based on whether the user wants the dense point cloud returned
  if(return.dense == FALSE){

    return(segmented)

  } else {
    ## get the indices of nearest sparse points
    s2d_lookup <- RANN::nn2(data=data.frame(X=segmented@data$X,Y=segmented@data$Y,Z=segmented@data$Z),
                            query=data.frame(X=las@data$X,Y=las@data$Y,Z=las@data$Z),
                            k=1)
    ## grab the treeID of the sparse point
    las@data$treeID<-sparse@data[s2d_lookup[[1]]]$treeID
    las@data$nnd<-s2d_lookup[[2]]

    ## assign the point if within subsampled sparse point density, otherwise 0
    las@data[las@data$nnd>(subsample.graph * 1.25)]$treeID<-0
    las@data[las@data$Classification==2]$treeID<-0

    # plot(las, color="treeID", backend="pcv")

    ## add the treeID values as attribute the dense point cloud
    las <- add_lasattribute(las, (las@data$treeID), name = "treeID", desc = "tree id")

    return(las)
  }

}

