#'Phenotyping Flow Cytometry Assays
#'
#'flowTypeFilter uses a simple threshold, polygone filters,
#'Kmeans, flowMeans or flowClust to split every channel to multiple partitions.
#'These partitions are then combined to generate a set of
#'multi-dimensional phenotypes.
#'
#'@param Frame A flowFrame (after transformation) that is going to be phenotyped.
#'
#'@param PropMarkers A vector of the indexes or names of the markers to
#'partition to specify phenotypes.
#'If `NULL`, all markers in the frame will be used.
#'When using `Filters` Method,
#'`PropMarkers `needs to be set for only channels that use numeric thresholds.

#'
#'@param MFIMarkers A vector of the indexes or names of the markers for which
#'MFIs must be measured.
#'If `NULL`, no markers will be used.
#'Currently it won't be calculated when filters are given.
#'
#'@param Methods  A single string specifying the method to use to determine
#'thresholds for partitioning of markers. Values can be "kmeans", "flowMeans",
#'"flowClust", "Thresholds", or "Filters". If "Thresholds"is specified,
#'user-specified thresholds must be provided via the Thresholds parameter.
#'If 'Filters' is specified, thresholds or polygon filters must be provided
#'as a list. These can be calculated using flowDensity or
#'flowCore `polygoneGate()`
#'
#'@param MarkerNames A vector of names for the channels. If `NULL`, the names
#'in `Frame` will be used.
#'
#'@param MaxMarkersPerPop An integer specifying the maximum number of markers to
#' use to define populations (how "deep" to phenotype). This should be less than
#'  or equal to `PropMarkers`.
#'  If `NULL`, will default to the length of PropMarkers.
#'
#'@param PartitionsPerMarker An integer or vector of integers specifing the
#'number of partitions per marker. If a single integer, this number will be used
#' for all markers. If a vector, the numbers will be matched with `PropMarkers`
#'  in order.
#'
#'@param Thresholds When `Method=="Thresholds"`, a list of vectors specifying
#'per-channel thresholds. Each list item corresponds to one marker,
#'and contains the threshold(s) for that marker. If only one vector is provided
#'in the list, then those thresholds will be used for all markers.
#'Otherwise, the list must be of the same length as PropMarkers.
#'When `Method=="Filters"`. a list of vectors/matrices. Each element corresponds
#'to either one (when it's a vector) or two (when it's a matrix) marker(s).
#'Matrices are gate boundaries that are compatible with polygoneGate in
#'flowCore,or flowDensity@filter. for each filter matrix, column names must be
#'channels corresponding to 2-D gate, that match FCS channel names.
#'Note: if `Methods == 'thresholds' | Methods == 'Filters'}`,
#'then `Thresholds` must be specified. If not, it is ignored.
#'
#'@param MemLimit Memory limit in GB. flowTypeFilter will do a sanity check
#'before executing, and if the total size of counts plus MFI values for
#'all populations would exceed `MemLimit`, will not run.
#'
#'@param verbose Boolean variable. If TRUE, information about different
#'processing tasks will be printed into the standard output.
#'
#'@param cores The number of cores to use for mclapply function when
#'`Methods == 'Filters'`
#'
#'@param ... For mclapply, optional arguments.
#'
#'@return New`PhenoType`object contation
#'`CellFreqs`: Object of class `Numeric` containing the cell frequencies
#'measured for each phenotype. Phenotype names are assigned as labels.
#'
#'@examples
#'##Load the library
#'library(flowTypeFilter)
#'data(DLBCLExample)
#'
#'##These markers will be analyzed
#'PropMarkers <- 3:5
#'MFIMarkers <- PropMarkers
#'MarkerNames <- c('FS', 'SS','CD3','CD5','CD19')
#'
#'##Run flowTypeFilter
#'Res <- flowTypeFilter(DLBCLExample, PropMarkers, MFIMarkers, 'kmeans',
#'MarkerNames);
#'
#'MFIs=Res@MFIs;
#'Proportions=Res@CellFreqs;
#'Proportions <- Proportions / max(Proportions)
#'names(Proportions) <- unlist(lapply(Res@PhenoCodes,
#'                                    function(x){return(decodePhenotype(
#'                                    x,Res@MarkerNames[PropMarkers],
#'                                    Res@PartitionsPerMarker))}))
#'
#'##Select the 30 largest phenotypes
#'index=order(Proportions,decreasing=TRUE)[1:30]
#'bp=barplot(Proportions[index], axes=FALSE, names.arg=FALSE)
#'text(bp+0.2, par("usr")[3]+0.02, srt = 90, adj = 0,
#'labels = names(Proportions[index]), xpd = TRUE, cex=0.8)
#'axis(2);
#'axis(1, at=bp, labels=FALSE);
#'title(xlab='Phenotype Names', ylab='Cell Proportion')
#'
#'##These phenotype can be analyzed using a predictive model
#'(e.g., classification or regression)
#'}
#'
#' @export

flowTypeFilter <- function(Frame,
                           PropMarkers=NULL,
                           MFIMarkers=NULL,
                           Methods='Filters',
                           MarkerNames=NULL,
                           MaxMarkersPerPop=NULL,
                           PartitionsPerMarker=2,
                           Thresholds=NULL,
                           MemLimit=4,
                           cores = getOption("mc.cores", 2L),
                           verbose=TRUE,...)
{
  ##############################################################################################################################
  # Argument processing
  ##############################################################################################################################
  ##If list of markers are not supplied, use all of the available channels
  if(is.null(PropMarkers))
  {
    if (Methods=="Filters" & (any (unlist(lapply(Thresholds,class))!="matrix")|any (unlist(lapply(Thresholds,class))!="data.frame") ))
      stop("You need to provide thresholds' channel numbers via 'PropMarkers' when using 'Filters' Method.")
    PropMarkers=c(1:length(exprs(Frame)[1,]));
    ind <- which(unlist(lapply(Thresholds,class))=="data.frame"|unlist(lapply(Thresholds,class))=="matrix")
    PropMarkers <- PropMarkers[-ind]
  }
  if(is.null(MFIMarkers))
    MFIMarkers=vector();


  ### Remove the row of FCS file name from Thresholds
  # if (any (unlist(lapply(Thresholds, class)) == "character"))
  # {
  #   ind <- which(unlist(lapply(Thresholds, class)) == "character")
  #   Thresholds <- Thresholds[-ind]
  # }

  ### Threshold filter should be a matrix
  if (Methods=="Filters" & any (unlist(lapply(Thresholds,class))!="matrix") )
  {
    ind <- which(unlist(lapply(Thresholds,class))=="data.frame"|unlist(lapply(Thresholds,class))=="matrix")
    if (length(ind)==0 | length(ind)!=(length(Thresholds)-length(PropMarkers)))
    {

      stop("You must provide all channel numbers of simple thresholds,
            and make sure your gating filters are of class 'matrix' or 'data.frame'.")

    }else
    {
      for ( i1 in ind)
        Thresholds[[i1]] <- as.matrix( Thresholds[[i1]])
    }
  }

  ## TODO: typecheck methods parameter.
  if (is.null(MaxMarkersPerPop) & Methods=="Filters")
  {
    MaxMarkersPerPop <- length(Thresholds)

  }else if (is.null(MaxMarkersPerPop))
  {
    MaxMarkersPerPop <- length(PropMarkers)
  }

  ##If list of markers are names, convert them to indexes.
  if (FALSE %in% is.numeric(PropMarkers))
    PropMarkers <- unlist(lapply(1:length(PropMarkers), function(i){which(PropMarkers[i]==colnames(exprs(Frame)))}))
  if (FALSE %in% is.numeric(MFIMarkers)&&length(MFIMarkers) > 0)
    MFIMarkers <- unlist(lapply(1:length(MFIMarkers), function(i){which(MFIMarkers[i]==colnames(exprs(Frame)))}))

  ##Parse method:

  VALID_METHODS <- c('Thresholds', 'Filters','flowMeans', 'kmeans', 'flowClust')
  Methods <- sapply(Methods, function(x){sub('thresh', 'Thresh', x)}) #fix case issue with thresholds
  Methods <- sapply(Methods, function(x){sub('filters', 'Filters', x)}) #fix case issue with filters

  if(! all(sapply(Methods, function(x){x %in% VALID_METHODS})) )
    stop(paste('Invalid method specified. Methods must be one of:', paste(VALID_METHODS, collapse = ', ')))

  if (length(Methods)!=1)
    stop('Only one method may be specified')

  if(length(Thresholds) == 1 & Methods=="Filters" & length(PropMarkers)>1)
  {
    stop('When providing Filters, either thresholds and/or filters need to be provided for channels.')

  }

  if( Methods=='Thresholds' && (!is.list(Thresholds)) )
    stop('Thresholds must be provided as a list of vectors.')

  if(length(Thresholds) == 1)
  {
    if(length(unique(PartitionsPerMarker)) > 1)
      stop('When markers have different numbers of partitions, you must specify Thresholds on a per-marker basis.')
    if(length(Thresholds[[1]])!=PartitionsPerMarker[1]-1)
      stop('When a single vector is provided for Thresholds, it must contain exactly PartitionsPerMarker-1 Thresholds.')
    Thresholds <- rep(Thresholds, length(PropMarkers))
  }




  # if (length(PartitionsPerMarker)!=length(PropMarkers))
  #  stop('PartitionsPerMarker must either be specified once for all markers, or be of the same length as PropMarkers.')


  if(length(Thresholds)==0 && (Methods=="Filters" | Methods=="Thresholds"))
    stop('When Thresholds or Filters is specified as a method, You must provide Thresholds via the "Thresholds" argument.')
  if (is.null(Thresholds))
    Thresholds <- vector("list",length(PropMarkers))

  ##If PartitionsPerMarker is a single value, replicate for all PropMarkers:
  if (length(PartitionsPerMarker)==1)
    PartitionsPerMarker<-rep(PartitionsPerMarker, length(Thresholds));
  #Get marker names if not provided:
  if(is.null(MarkerNames))
    MarkerNames <- as.vector(Frame@parameters@data$name)



  ##############################################################################################################################
  # Memory check
  ##############################################################################################################################

  ##Before even starting, perform a sanity check for memory usage
  NumPops <- calcNumPops(PartitionsPerMarker, MaxMarkersPerPop)

  MemUse <- calcMemUse(NumPops, length(Thresholds), length(MFIMarkers), nrow(Frame), MaxMarkersPerPop, max(PartitionsPerMarker)) / 10^9
  if (verbose)
    message(sprintf('Estimated memory required: %f GB of RAM', MemUse))
  if(MemUse > MemLimit)
    stop(paste('Calling flowType with these parameters would require', MemUse, 'GB of RAM, but MemLimit is', MemLimit, 'GB.\n Try reducing MaxMarkersPerPop or the number of MFIMarkers.'))

  ##############################################################################################################################
  # Perform partititoning
  ##############################################################################################################################


  ##express the flowFrame and get channel numbers
  if (Methods=="Filters")
  {
    names(Thresholds)[which(unlist(lapply(Thresholds,class))!="matrix")] <- colnames(exprs(Frame))[PropMarkers]
    names(Thresholds)[which(unlist(lapply(Thresholds,class))=="matrix")] <-"gate"

    ##Fix the joint mark difference of the Threshold names
    inds <-which(names(Thresholds)=="gate")
    for (i in inds) {
      colnames(Thresholds[[i]]) <- gsub(".", "-",colnames(Thresholds[[i]]), fixed = TRUE)
    }

    ## Rearrange the Thresholds and make gates at the end
    gates <- Thresholds[inds]
    Thresholds <- append(Thresholds[-inds], gates)

    ## Get all channels related to thresholds and express the Frame
    gate.channels <- find.channels(Thresholds, Frame)
    all.channels <- c(PropMarkers, gate.channels)
    X <- exprs(Frame)[, sort(unique(all.channels))]

    ## Get the number lables of markers correspoding to X channels, passing them to C code
    ThresChannels <- match(all.channels, sort(unique(all.channels)))
  }else{
    X <- exprs(Frame)[,PropMarkers]
  }
  ##number of markers
  M <- length(Thresholds)
  ##number of events
  N <- nrow(X)
  ##the matrix that will hold the partitions
  Partitions <- matrix(0,M,N);

  ##Use the respective method.

  #K-means partitioning:
  if (Methods=='kmeans'){
    for (i in 1:M){
      km<-kmeans(X[,i], PartitionsPerMarker[i], nstart=50)$cluster

      ##Sort clusters in increasing order by mean:
      means <- sapply(unique(km), function(x){mean(X[which(km==x),i])})
      names(means) <- unique(km)
      means <- sort(means)
      new.km <- rep(0, length(km))
      for (cluster.ind in 1:PartitionsPerMarker[i])
      {
        to.replace <- which(names(means)==cluster.ind)
        new.km[which(km==cluster.ind)] <- to.replace
      }
      Partitions[i,] <- new.km;
    }
  }

  #flowMeans partitioning:
  if (Methods=='flowMeans'){
    for (i in 1:M){
      km<-flowMeans(X[,i], NumC=PartitionsPerMarker[i], MaxN=10, nstart=10)@Label
      ##Sort clusters in increasing order by mean:
      means <- sapply(unique(km), function(x){mean(X[which(km==x),i])})
      names(means) <- unique(km)
      means <- sort(means)
      new.km <- rep(0, length(km))
      for (cluster.ind in 1:PartitionsPerMarker[i])
      {
        to.replace <- which(names(means)==cluster.ind)
        new.km[which(km==cluster.ind)] <- to.replace
      }

      Partitions[i,] <- new.km;
    }
  }

  if (Methods=='flowClust'){
    for (i in 1:M){
      res <- flowClust(Frame, varNames=colnames(exprs(Frame))[i], K = PartitionsPerMarker[i], level=1);

      km=map(res@z);
      km=replace(km, which(is.na(km)), 1)

      ##Sort clusters in increasing order by mean:
      means <- sapply(unique(km), function(x){mean(X[which(km==x),i])})
      names(means) <- unique(km)
      means <- sort(means)
      new.km <- rep(0, length(km))
      for (cluster.ind in 1:PartitionsPerMarker[i])
      {
        to.replace <- which(names(means)==cluster.ind)
        new.km[which(km==cluster.ind)] <- to.replace
      }

      Partitions[i,] <- new.km;
    }
  }

  #If method is Thresholds, still calculate partition membership for later plotting and labelling purposes:
  if (Methods=='Thresholds'){
    for (Marker in 1:length(Thresholds))
    {
      #Calculate partition membership for each partition for this marker:
      marker.vec <- rep(1, ncol(Partitions)) #Set to 1 (first partition) by default
      for(partition in 2:(length(Thresholds[[Marker]])+1))
      {
        #Then go through each partition, setting everything greater than its threshold
        #Partition membership gets progressively overwritten until all are right
        marker.vec[which(X[,Marker] >= Thresholds[[Marker]][partition-1])] <- partition
      }
      Partitions[Marker,] <- marker.vec
    }
  }


  if (Methods=='Filters'){
    for (Marker in 1:length(Thresholds))
    {

      if(names(Thresholds)[Marker]=="gate")
      {
        #Calculate partition membership for each partition for this marker:
        marker.vec <- rep(1, ncol(Partitions))
        channels <- colnames(Thresholds[[Marker]])
        pts <- sp::point.in.polygon(point.x = X[,channels[1]], point.y = X[,channels[2]],
                                    pol.x = Thresholds[[Marker]][,1],pol.y =Thresholds[[Marker]][,2])
        #Then go through each partition, setting everything greater than its threshold
        #Partition membership gets progressively overwritten until all are right
        marker.vec[c(which(pts==1),which(pts==2),which(pts==3))] <- 2
        Partitions[Marker,] <- marker.vec
      }else{
        marker.vec <- rep(1, ncol(Partitions)) #Set to 1 (first partition) by default
        channel <- names(Thresholds)[Marker]
        for(partition in 2:(length(Thresholds[[Marker]])+1))
        {
          #Then go through each partition, setting everything greater than its threshold
          #Partition membership gets progressively overwritten until all are right
          marker.vec[which(X[,channel] > Thresholds[[Marker]][partition-1])] <- partition
        }
        Partitions[Marker,] <- marker.vec
      }
    }
  }


  ###Turn partitions into Thresholds to pass down to CPP code:
  partToThresh <- function(ThisChan, Partitions, PropMarkers, ThisExpr)
  {
    ThisPart <- t(Partitions)[, ThisChan]
    PartLabels <- unique(ThisPart)

    #Find maxima of all partitions:
    MaxParts <- sapply(PartLabels, function(x){max(ThisExpr[which(ThisPart==PartLabels[x]),ThisChan])})
    Thresholds <- sort(MaxParts)[1:(length(PartLabels)-1)]

    return(Thresholds)
  }

  if(Methods %in% c('flowMeans', 'kmeans', 'flowClust'))

    Thresholds <- lapply(1:length(PropMarkers), partToThresh, Partitions, PropMarkers, X)
  for (i in 1:M){
    ##because flowMeans's 1D Thresholds are not perfectly aligning, we calculate the partitions one more time after the Thresholds are finalized.
    if (Methods=='flowMeans'){
      for (Marker in 1:length(Thresholds))
      {
        marker.vec <- rep(1, ncol(Partitions))
        for(partition in 2:(length(Thresholds[[Marker]])+1))
        {
          marker.vec[which(X[,Marker] >= Thresholds[[Marker]][partition-1])] <- partition
        }
        Partitions[Marker,] <- marker.vec
      }
    }
  }



  ##############################################################################################################################
  # Calculate MFIs and counts via CPP code and return
  ##############################################################################################################################

  #Call down  to cpp code to calculate:
  if(length(MFIMarkers) > 0)
  {
    MFIData <- matrix(exprs(Frame)[,MFIMarkers], ncol=length(MFIMarkers))
  }else
  {
    MFIData <- matrix()
  }
  if (Methods=="Filters"){
    isFilter <- TRUE
    res <- .Call('countCells', isFilter, as.integer(PartitionsPerMarker), Thresholds, as.integer(MaxMarkersPerPop), ThresChannels, MFIData, X, NumPops, verbose)
  }else{
    isFilter <- FALSE
    res <- .Call('countCells', isFilter, as.integer(PartitionsPerMarker), Thresholds, as.integer(MaxMarkersPerPop), PropMarkers, MFIData, X, NumPops, verbose)
    names(Thresholds) <- colnames(X)
  }
  print("Returning counts...")
  names(res$counts) <- names(res$codes) <- NULL
  Counts <- res$counts;
  if (length(MFIMarkers)>0){
    Means <- res$mfis;
    Means[which(Counts==0),]=NA
  }else
  {
    Means=matrix();
  }
  Codes <- res$codes;

  if (length(MFIMarkers)>0){

    colnames(Means)=colnames(exprs(Frame))[MFIMarkers];
  }
  #If still NULL, just make it be channel names of the Frame:
  if(is.null(MarkerNames))
    MarkerNames <- colnames(X)

  Partitions=t(Partitions);
  return(new("Phenotypes", CellFreqs=Counts, PhenoCodes=Codes, MFIs=Means, PropMarkers=PropMarkers, MFIMarkers=MFIMarkers, MarkerNames=MarkerNames, Partitions=Partitions, PartitionsPerMarker=PartitionsPerMarker, Thresholds=Thresholds));
}



