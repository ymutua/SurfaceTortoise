#' @title SurfaceTortoise

#' @description Optimizing spatial sampling using the Surface Tortoise
#' algoritm. Grid sampling and random  sampling are also available. All three sampling
#' designs can optionally be stratified by a square grid to ensure spatial coverage.
#'
#' @author Kristin Piikki, Mats Söderström & John Mutua ,  \email{kristin.piikki@@slu.se}

#' @param x Raster dataset. Required for method = directed. The raster must have
#' a defined coordinate system and must be of class numeric. If x is a raster stack or
#' raster brick, the first principal component of the multiple layers will be
#' used for sampling. If the raster dataset has a single layer, it will be used
#' as is.

#' @param y SpatialPolygonsDataframe delineating the area to be sampled. Required
#' for method = 'grid' and method = 'random'. Optional for method = 'directed.
#' The SpatialPolygonsDataframemust must have a defined coordinate system and,
#' if a raster is provided, the coordinate system shall be the same as for the
#' raster. If x and y are not completely overlapping, their intersection will be
#' sampled.

#' @param method Sampling method: 'directed' = directed sampling (SurfaceTortoise
#' algorithm), 'grid' = regular sampling (center points of strata) and
#' 'random' = random points. Default is 'directed'
#'
#' @param filter An integer. Side of the square window (number of
#' raster cells, original resolution) used for mean filtering of the raster.
#' Default = 1 (no filtering)
#'
#' @param resolution A number. If provided, the raster data vill be resampled
#' to this resolution. Optional.
#'
#' @param p_idw An integer. Power exponent used for idw-interpolation
#' (method = 'directed'). Default is 2.
#'
#' @param nmax_idw An integer. Number of neighbouring samples used for
#' idw-interpolation (method = 'directed'). Default is 8.
#'
#' @param edge A number. Buffer zone (metre) inside the sampled area border,
#' where sampling is prohibited. Optional.
#'
#' @param strat_size A number. Cell side (metre) of a square stratification grid.
#' Optional. #' If both strat_size and stop_n are specified. stop_n overruns
#' this argument #' with an adjusted strat_size. If strat_size is not specified.
#' The sampling will be done without stratification. If strat_size = 0,
#' stratification size will be comuted from the number of samples. Negative
#' values are not allowed.
#'
#' @param min_dist A positive number. Minimum distance allowed between samples.
#' Valid for the 'random' and the 'directed' methods.
#'
#' @param stop_n An integer. The number of samples to place. If not provided,
#' it will be conuted from the numbers of strata generated from the specificed
#' stratication size (argument strat_size) and the number of samples to place
#' per stratum (argument stop_dens).
#'
#' @param stop_dens An integer. The number of samples to place in each stratum.
#' Does not apply for method = 'grid' (always stop_dens = 1) and not for
#' non-stratified sampling. Default is 1.
#'
#' @param plot_results Logical. Shall results be plotted? Default is FALSE.

#' @return  A list with
#' 1) sampled_raster = the sampled raster (only if method = 'directed')
#' 2) samples = a spatialPointsDataFrame with sample locations
#' 3) sampled_area = a SpatialPolygonsDataFrame with a polygon for the sampled area.
#' 4) stratification = a a SpatialPolygonsDataFrame with the stratification polygons.
#' 5) feedback= a dataframe with generated text messages.

#' @details The Surface Tortoise algorithm for directed sampling
#' uses a raster dataset to find optimal sample locations.
#' The sampling strategy is based on the principle that an interpolation of the
#' samples should be as similar as possible to the guide raster. When sample
#' locations are identified, first the center point of the raster cell with the
#' maximum deviation from the covariate raster mean is sampled. Then the raster
#' cell with the maximum deviation from the first sampled raster cell is sampled.
#' From then on, the values of the sampled raster cells are interpolated by
#' inverse distance weighting (idw) and the center point of the raster cell with
#' the largest absolute difference to the guide raster (error) is sampled.
#' A new idw interpolation is made and a new cell is sampled. This is repeated
#' is reached.The sampling can be stratified by a square grid. When a sample has
#' been placed in a stratum, no more samples will be placed in that stratum again
#' until all other strata have been sampled. The likelihood for a clipped stratum,
#' e.g. at the edge of the area to be sampled, is equal to the area of that stratum
#' divided by the area of a full stratum.
#'
#' The optional raster processing steps:
#' (is done) is carried out in the folowing order:
#' 1) mean filtering (argument: filter)
#' 2) resampling to specified resolution (argument: resolution),
#' 3) computation of first pricipal component (if x is a rastr stack or raster
#' brick with multiple layers).

#' @import raster
#' @import gstat
#' @import rgeos
#' @import sp

#' @references Olsson, D. 2002. A method to optimize soil sampling from
#' ancillary data. Poster presenterad at: NJF seminar no. 336,
#' Implementation of Precision Farming in Practical Agriculture, 10-12 June
#' 2002, Skara, Sweden.

#' @examples
#' data(boundary)
#' grid.sampling<-tortoise(y=boundary,method='grid',edge=30,strat_size=50,
#' min_dist=10,plot_results=FALSE)

#' @export
tortoise<-function(x=NULL,
                   y=NULL,
                   method='directed',
                   edge=0,
                   strat_size=NULL,
                   min_dist=0,
                   p_idw=2,
                   nmax_idw=8,
                   resolution= NULL,
                   filter=1,
                   stop_n=NULL,
                   stop_dens=1,
                   plot_results=F
                   )
  {
  ###FEEDBACK FUNCTION
  ###FEEDBACK FUNCTION
  ###FEEDBACK FUNCTION

  feedback<-'Messages'
  fback<-function(t){print(t); return (c(feedback,t))}


  ###PREPARE DATA
  ###PREPARE DATA
  ###PREPARE DATA

  #are data provided?
  if(method=='directed'){
    if(is.null(x)) {
      stop(call. =F, 'Raster data not provided.')
    }

  }
  if(method!='directed'){
    if(is.null(y)) {
      stop(call. =F, 'Boundary polygon not provided.')
    }
  }
  #are data projected?
  if(method=='directed')if(is.na(crs(x))){
    stop(call. =F, 'Raster data not projected.')
  }
  if(!is.null(y)) if(is.na(crs(y))) {
    stop(call. =F, 'Boundary polygon not projected.')
    }

  #are the raster and polygon data projected onto the same coordinate system?
  if(method=='directed' & !is.null(y)){
    if(compareCRS(x, y)==F) {
      t<-'The coordinate systems of the raster data and the polygon data are not the same. The polygon data will be projected to the coordinate system of the raster data.'
      feedback<-fback(t)
      y<-spTransform(x=y, CRSobj=crs(x))
    }
  }

  #is the polygon data in a polar or cartesian coordinate system? Reproject to
  #web mercator if the latter
  if(!is.null(y))if (is.projected(y)==F) {
    crs_y<-crs(y)
    web_mercator<-'+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs'
    y<-spTransform(x=y, CRSobj=CRS(web_mercator))
    if(method == 'directed')x<-spTransform(x=x, CRSobj=CRS(web_mercator))
    t<-'The  data were reprojected to the wgs84 web mercator (auxiliary sphere) coordinate system.'
    feedback<-fback(t)
  }

  #are the rasters of class numeric?
  if(method == 'directed') if(!is.numeric(values(x))){
    t<-'Error: the raster(s) is/are not of class numeric.'
    feedback<-fback(t)
  }

  #is the polygon overlapping the raster?
  if(!is.null(x) & !is.null(y) & (method=='directed')){
    if(!ncell(crop(x,y))>0) stop(call. =F, 'x and y are not overlapping')
  }

  #is a too high number of samples per strata (stop_dens) specified?
  if(stop_dens>10){
    t<-'Warning: A large number of samples per strata has been specified. The sampling may take a very long time. Consider changing argument stop_dens to a smaller value or use the default.'
    feedback<-fback(t)
  }
  #set stop_dens to 1 if method = grid
  if (method=='grid')stop_dens<-1

  #is a too large total number of samples (stop_n) specified
  if(!is.null(stop_n))if(stop_n>500) {
    t<-'Warning: A large number of samples has been specified. The sampling may take a very long time. Consider changing argument stop_n to a smaller value or use the default.'
    feedback<-fback(t)
  }

  #is strat_size == 0 when stop_n has not been specified
  if(is.null(strat_size)& is.null(stop_n))  {
    stop(call. =F, 'At least one of the arguments stop_n and strat_size needs to be specified')
  }
  if(!is.null(strat_size)) if(strat_size ==0 & is.null(stop_n))  {
    stop(call. =F, 'If number of samples is not specified (argument: stop_n), strat_size must be >=1')
  }

  #store projection
  prj<-crs(x)

  #create polygon(y), if not provided
  if(!is.null(x) & is.null(y)){
    bbx<-rasterToPolygons(x, dissolve=TRUE)
    crs(bbx)<-crs(x)
    y<-bbx
  }

  #clip the polygon, if the raster does not fill it
  if(!is.null(x) & method == 'directed'){
    bbx<-rasterToPolygons(x, dissolve=TRUE)
    crs(bbx)<-crs(x)
    y<-rgeos::gIntersection(y, bbx, byid=F)

  }

  #shrink polygon, if a buffer (argument edge) is specified
  y_large<-y #save a copy

  if(edge>0) {if (!is.null(buffer(y,-edge))) {
    y <- buffer(y,-edge)} else {
      t<-'Error: The buffer is too large, please chose a smaller value of argument edge.'
      feedback<-fback(t)
    }
  }

  #aggregate features, if y contains multiple features (rows)
  if(!is.null(nrow(y))) if(nrow(y)>1) y <- aggregate(y,dissolve=T)

  #roughly crop and mask raster by the polygon
  if(!is.null(x) & method!='grid'& method!='stratrand'){
    x<-crop(x, buffer(y, 5*res(x)[1]))
  }

  #mean filter raster data
  if(method=='directed' & filter>1) {
    aa<-is.na(x)
    if(filter %% 2 == 0) filter<-filter+1 #if filter is even add one
    w<-matrix(1, filter, filter)
    if(nlayers(x)==1){
      x<-raster::focal(x,w=w, fun=mean, na.rm=TRUE, NAonly = TRUE, pad=T)
    }
    if(nlayers(x)>1){
      for (i in 1:nlayers(x)){
        x[[i]]<-raster::focal(x[[i]],w=w, fun=mean, na.rm=TRUE, NAonly = TRUE, pad=T)
      }
    }
    x[aa]<-NA
    t<-'Your raster have been mean filtered. See argument filter.'
    feedback<-fback(t)
  }

  #resample if raster resolution is too fine and no new resolution is specified
  finest_resolution<- signif((sqrt(area(y)))/300, 1)
  if(method == 'directed' & is.null(resolution)) if(res(x)[1] < finest_resolution){
    template<-raster(res=finest_resolution, ext=extent(x), crs=prj)
    x<-resample (x=x, y=template, method='bilinear')
    t<-paste0('Your raster was resampled to a cell size of ', finest_resolution, ' m to speed up the processing. If you do not want this to happen, please provide the desired resolution (using argument resolution).')
    feedback<-fback(t)
  }
  #warn if specified resolution is too fine
  if(method == 'directed' & !is.null(resolution))if(resolution < finest_resolution){
    t<-paste0('The specified resolution was very fine in relation to the size of the area. Processing will take a very long time. Consider setting a coarser resolution to speed up the processing.')
    feedback<-fback(t)
  }
  #resample to the resolution specified by the user
  if (method == 'directed' & !is.null(resolution)) if(res(x)[1]!=resolution) {
    template<-raster(res=resolution, ext=extent(x), crs=prj)
    x<-resample (x=x, y=template, method='bilinear')
    t<-paste0('Your raster was resampled to the specified resolution of ', resolution, ' m .')
    feedback<-fback(t)
  }
  #is raster resolution too coarse in relation to other sampling specifications
  if(method=='directed') if (stop_n>area(y)/((res(x)[1])^2)){
    t<-'The raster resolution is very coarse in relation to the size of the area and the number of samples to place. Consider increasing the area to be sampled, reduce the numbers of samples to place (argument: stop_n) or resample the raster to another resolution (argument: resolution).'
    feedback<-fback(t)
  }

  #final cropping and masking of raster(s) by polygon
  if(!is.null(x) & method!='grid'& method!='stratrand'){
    x<-crop(x, y)
    x<-mask(x=x, mask=y)
  }

  #compute pca of rasterstack
  if(method=='directed') if(nlayers(x)>1){
    df <- data.frame(values(x))
    cc<-stats::complete.cases(df)
    pc1 <- data.frame(stats::prcomp(df[cc,], scale = TRUE)['x'])[,1]
    df<- data.frame(coordinates(x)[cc,], pc1)
    names(df) <- c ('x', 'y', 'pc1')
    coordinates(df)<-~x + y
    x<-rasterize(df, x)[[2]]
    crs(x)<-crs(y)
    t<- 'More than one raster layer was provided so the first pricipal component of them will be used to direct the sampling.'
    feedback<-fback(t)
  }

  ###CREATE STRATIFICATION
  ###CREATE STRATIFICATION
  ###CREATE STRATIFICATION

  #if strat_size is not provided, find the approximate strat_size fulfilling stop_dens and stop_n
  if(is.null(strat_size) & method == 'grid') strat_size<-0
  if(!is.null(strat_size))if(strat_size==0 | !is.null(stop_n)) {
    autostratification<-T
    if (!is.null(stop_n)){
      t<-'Argument strat_size was ignored. The size of the stratification grid was computed from the the total number of samples (argument: stop_n) and the numbers of samples per stratum (argument: stop_dens).'
      feedback<-fback(t)
    }
    }else autostratification<-F
  if(is.null(strat_size)) autostratification<-F
  if(autostratification) strat_size<-signif(sqrt(area(y))/sqrt(stop_n/stop_dens), 2)

  #is a too fine stratification grid used?
  if(!is.null(strat_size) & autostratification==F) {
    check<-area(y)/(strat_size^2)>500
    if(check){
      t<-'Warning: A very fine stratification grid is specified in relation to the size of the area. The sampling will take very long time. Consider changing argument strat_size to a larger value or use the default.'
      feedback<-fback(t)
    }
  if (method =='directed'){
    check<-strat_size<res(x)[1]
    if(check){
      t<-'Warning: The stratification grid has a smaller cell size than the raster. Please change argument strat_size to a larger value.'
      feedback<-fback(t)
    }
  }
  }

  #create a stratification grid
  if (!is.null(strat_size)){
    #define funtion for stratification
    f<-function(p, s, d){
      e<-extent(p)+c(-s, s,-s, s )
      strat<-raster(ext = e, resolution = s, crs= crs(p))
      values(strat)<-1:ncell(strat)
      strat<-rasterToPolygons(strat)
      names(strat)<-'id'
      if(method != 'grid') strat<-intersect(strat, p)
      strat$n<-0
      strat$target_n<-round(d*(area(strat)/max(area (strat))))
      strat$round_remain_n<-0
      return(strat)
    }
    strat<-f(p=y, s=strat_size, d=stop_dens)

    #if grid sampling, iteratively adjust strat_size to reach stop_n
    if(method == 'grid' & !is.null(stop_n)){
      for (x in 1:100){
        if(x==1)n<-sum(strat$target_n)
        if(stop_n>n)strat_size<-strat_size*0.9995
        if(stop_n<n)strat_size<-strat_size/0.9995
        if(stop_n==n)break
        strat<-f(p=y, s=strat_size, d=1)
        p.sp<-gCentroid(strat, byid=T)
        p.sp<-p.sp[y,]; n<-length(p.sp@coords[,1]) #ugly workaround to fin number of samples
        }
    if(stop_n!=n){
      t<-'Warning: Could not find optimal stratificaiton size by iteration. the number is samples will not be as specified by agrument stop_n. Please try to specify an approximate stratification size (argument: strat_size) to find a solution.'
      feedback<-fback(t)
    }
    t<-paste0('The stratification grid has been created. The cell size is ', round(strat_size), ' m * ', round(strat_size),' m.')
    feedback<-fback(t)
    }


    }
  #compute stop_n if not specified
  if(is.null(stop_n)){
      stop_n<-sum(strat$target_n)
      t<-paste0('The number of samples to place has been determined automatically.')
      feedback<-fback(t)
    }
    if(autostratification & sum(strat$target_n) < stop_n){
      strat$target_n<-round((stop_dens+1)*(area(strat)/max(area (strat))))
    }


  #if a non-stratified method is used, just use the (buffered, polygon) as strat (one stratum)
  if (is.null(strat_size)) {
    strat<-y
    strat$id<-1
    strat$n<-0
    strat$target_n<-stop_n
    strat$round_remain_n<-0
    t<-'A bounding polygon for the sampling has been created.'
    feedback<-fback(t)
  }

  ###DIRECTED SAMPLING
  ###DIRECTED SAMPLING
  ###DIRECTED SAMPLING

  if (method == 'directed'){

    #identify number of samples to take
    if(is.null(stop_n)){
      stop_n<-sum(strat$target_n)
      t<-paste0('The number of samples to place has been computed. It is ', stop_n)
      feedback<-fback(t)
    }

    #run sampling loop
    for (i in 1:stop_n){

      #computes samples left to take in strata in total
      strat$tot_remain_n <-strat$target_n-strat$n

      #refill no of samples left to take in this round, if needed
      if(sum (strat$round_remain_n)==0) strat[strat$tot_remain_n>0, 'round_remain_n']<-1

      #calculate prediction raster (pred)
      if(i==1) {
        pred<- x
        pred[!is.na(x)]<-mean(x@data@values, na.rm=T)
      }
      if (i>1) {
        mod <- gstat(id = 'r', formula = r~1, locations = p.sp, nmax=nmax_idw, set=list(idp=p_idw))
        pred <- interpolate(x, mod, debug.level = 0)
        pred <- mask(pred, x)
      }

      #calculate aboslute error raster (ae)
      ae<-abs(x-pred)

      #save mae in p.sp
      if(i>1){p.sp[(i-1), 'mae']<-mean(ae@data@values, na.rm=T)}

      #mask sampling area 1: keep only strata with samples left to take
      m<-strat[strat$round_remain_n==1,]
      if(nrow(m)>0) ae<-mask(ae, m)

      #mask sampling area 2: exclude areas wihtin min_distance to points
      if (nrow(m)>0& i>1 & min_dist > 0) ae<-mask(ae, buffer (p.sp, min_dist), inverse=T)    #test if there is any area left to sample

      #test if there is any area left to sample
      if(cellStats(!is.na(ae), stat='sum')==0) {
        t<-'Error: The specified number of samples could not be placed. No area left to sample. Try a smaller min_dist, a smaller edge and/or a smaller stop_n!'
        feedback<-fback(t)
        break
      }


      #identify cell with maximum absolute error (ae), or for the first sample (i=1) the minimum ae
      if(i==1) {
        cell<-which.min(ae)
        if (length(cell)>1) cell<-sample(x=cell, size=1)
      }
      if(i>1) {
        cell<-which.max(ae)[1]
        if (length(cell)>1) cell<-sample(x=cell, size=1)
      }
      pi.r<-rasterFromCells(ae, cell, values=T)
      pi.sp<-rasterToPoints(pi.r, spatial=T); crs(pi.sp)<-crs(x)
      pi.sp$n<-i
      pi.sp$r<-raster::extract(x=x, y=pi.sp)
      pi.sp$mae<-NA  #create column
      if(i==1)p.sp<-pi.sp
      if(i>1)p.sp<-rbind(p.sp, pi.sp)

      #update strat
      id<-as.vector(strat[pi.sp, 'id']@data)[1,1]
      d<-data.frame(strat@data)
      d[d$id==id, 'n']<-d[d$id==id, 'n']+1
      d[d$id==id, 'round_remain_n']<-d[d$id==id, 'round_remain_n']-1
      strat@data<-d

      #plot results
      if(plot_results==T){
        graphics::par(mfrow=c(2,2))
        plot(y_large, main= 'Original surface')
        plot(x, legend=F, add=T,  )
        plot(strat, add=T)

        plot(y_large, main= 'Reconstructed surface')
        plot(pred, legend=F, add=T)
        plot(strat, add=T)
        plot(p.sp, add=T, pch=20)

        plot(y_large, main= 'Difference between the surfaces')
        plot((pred-x), add=T,legend=F)
        plot(strat, add=T)
        plot(p.sp, add=T, pch=20)

        if(i>3) {
          ymx<-p.sp@data[3,'mae']*2

          plot(p.sp$n[3:nrow(p.sp)], p.sp$mae[3:nrow(p.sp)],
               ylab='Difference', xlab='No of samples',
               xlim=c(0, stop_n),  ylim=c(0, ymx),
               main= 'Difference vs no of samples'
          )
        }
      }
    }

    #give feedback
    t<-'Sampling ready.'
    feedback<-fback(t)
  }

  ###GRID SAMPLING
  ###GRID SAMPLING
  ###GRID SAMPLING

  if (method == 'grid'){

    #find centroids of strata
    p.sp<-gCentroid(strat, byid=T) ##create sentroids to be used for grid sampling
    crs(p.sp)<-crs(y)


    if(class(p.sp)=="SpatialPoints"){
      p.sp.df <- data.frame(a1 = rep(1:5, length.out=length(p.sp)))
      p.sp <- SpatialPointsDataFrame(p.sp, p.sp.df)
    }

    p.sp<-intersect(p.sp, y)
    strat<-intersect(strat, y)
    p.sp$id<-1:nrow(coordinates(p.sp))

    #plot results
    if(plot_results==T){
      graphics::par(mfrow=c(1,1))
      plot(y_large,
           main='Grid sampling')
      plot(strat, add=T)
      plot(p.sp, pch=20, add=T)
    }

    #rename columns
    p.sp<-p.sp[,'id']
    names(p.sp)<-'n'

    #give feedback
    t<-'Sampling completed.'
    feedback<-fback(t)
  }

  ###RANDOM SAMPLING
  ###RANDOM SAMPLING
  ###RANDOM SAMPLING

  if (method == 'random'){

    #identify number of samples to take
    if(is.null(stop_n)){
      stop_n<-sum(strat$target_n)
      t<-paste0('The number of samples to place has been computed. It is ', stop_n)
      feedback<-fback(t)
    }

    #run sampling loop
    for (i in 1:stop_n){

      #computes samples left to take in strata in total
      strat$tot_remain_n <-strat$target_n-strat$n

      #refill no of samples left to take in this round, if needed
      if(sum (strat$round_remain_n)==0) strat[strat$tot_remain_n>0, 'round_remain_n']<-1

      #mask sampling area 1: keep only strata with samples left to take
      m<-strat[strat$round_remain_n==1,]

      #mask sampling area 2: exclude areas wihtin min_distance to points
      if(nrow(m)>0 & i>1 & min_dist > 0) {
        cl<-buffer (p.sp, min_dist)
        temp <- gDifference(m, cl, byid = TRUE)
        if(!is.null(temp)){
          row.names(temp) <- row.names(m@data)  #workaround
          m <- SpatialPolygonsDataFrame(temp, data = m@data)
        }else{
          t<-'Error: The specified number of samples could not be placed. No area left to sample. Try a smaller min_dist, a smaller edge and/or a smaller stop_n!'
          feedback<-fback(t)
          break
        }
      }


      #place a random sample
      pi.sp<-sp::spsample(x=m, n =1, type='random'); crs(pi.sp)<-crs(y)
      pi.sp$n<-i
      if(i==1)p.sp<-pi.sp else p.sp<-rbind(p.sp, pi.sp)

      #update strat
      id<-as.vector(strat[pi.sp, 'id']@data)[1,1]
      d<-data.frame(strat@data)
      d[d$id==id, 'n']<-d[d$id==id, 'n']+1
      d[d$id==id, 'round_remain_n']<-d[d$id==id, 'round_remain_n']-1
      strat@data<-d

      #plot results
      if(plot_results==T){
        graphics::par(mfrow=c(1,1))
        plot(y_large, main='Random stratified sampling')
        plot(strat, add=T)
        plot(p.sp, pch=20, add=T)
      }
    }

    #give feedback
    t<-'Sampling completed.'
    feedback<-fback(t)
  }

  ###COMPILE RESULTS
  ###COMPILE RESULTS
  ###COMPILE RESULTS

  #prepare  data
  if (method=='directed') sampled_raster<-x
  sampled_area<-y
  stratification<-strat[,c('id', 'n')]
  if (method=='directed') samples<-p.sp[,c('n', 'mae')] else samples<-p.sp[,'n']

  #list all results
  if (method=='directed'){
    results<-list(sampled_raster, sampled_area, stratification, samples, feedback)
  } else {
    results<-list(sampled_area, stratification, samples, feedback)
  }

  #return
  return(results)
  }

