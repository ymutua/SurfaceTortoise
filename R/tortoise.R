#' @title tortoise

#' @description Optimizing spatial sampling using the Surface Tortoise
#' algorithm. Grid sampling and random  sampling are also available. All three sampling
#' designs can optionally be stratified by a square grid to ensure spatial coverage.

#' @param x SpatRaster. Required for method = directed. The raster must have
#' a defined coordinate system, which is projected, and must be of class numeric.
#' If x has more than onelayer, the first principal component of the layers will be
#' used for sampling.

#' @param y SpatVector delineating the area to be sampled. Required
#' for method = 'grid' and method = 'random'. Optional for method = 'directed.
#' The SpatVector must have a defined coordinate system, which is projected.
#' If a SpatRaster is #' provided, the coordinate system shall be the same as
#' for the raster. If x and y are not completely overlapping, their intersection will be
#' sampled.

#' @param method Sampling method: 'directed' = directed sampling (SurfaceTortoise
#' algorithm), 'grid' = regular sampling (center points of strata) and
#' 'random' = random points. Default is 'directed'.
#'
#' @param filter A positive integer. Side of the square window (number of
#' raster cells, original resolution) used for mean filtering of the raster.
#' Default = 1 (no filtering).
#'
#' @param resolution A positive number. Optional. If provided, the raster data
#' will be resampled to this resolution.
#'
#' @param p_idw An integer. Exponent used for idw-interpolation
#' (method =='directed'). Default = 2.
#'
#' @param nmax_idw An integer. Number of neighbouring samples used for
#' idw-interpolation (method =='directed'). Default = 8.
#'
#' @param edge A positive number. Optional. Buffer zone (unit of the coordinate
#' reference system) inside the sampled area border,
#' where sampling is prohibited.
#'
#' @param strat_size A positive number, 0, or NULL. Optional. Cell side  of a square
#' stratification grid (unit of the coordinate reference system).If both
#' strat_size and stop_n are specified, stop_n overruns this argument with an
#' adjusted strat_size. If strat_size is NULL or 0, the sampling will be done
#' without stratification.
#'
#' @param min_dist A positive number. Optional. Minimum distance allowed between
#' samples. Valid for the 'random' and the 'directed' methods. one can for example
#' set min_dist to the diameter of the sample support to prevent overlapping
#' samples.
#'
#' @param stop_n A positive integer, or NULL. Optional. The number of samples to place.
#' If NULL, it will be computed from the numbers of strata generated from the
#' specified stratification size (argument strat_size) and the number of samples
#' to place per stratum (argument stop_dens).
#'
#' @param stop_dens A positive integer. The number of samples to place in each stratum.
#' Does not apply for method = 'grid' (always stop_dens = 1) and not for
#' non-stratified sampling. Default = 1.
#'
#' @param plot_results Logical. Shall results be plotted? Default is FALSE.

#' @return  A list with
#' 1) sampled_raster = the sampled raster (only if method =='directed')
#' 2) samples = a SpatVector of points (the sample locations)
#' 3) sampled_area = a SpatVector of polygons (the sampled area).
#' 4) stratification = a SpatVector of polygons (stratification grid).
#' 5) feedback = a dataframe with generated text messages.

#' @details The Surface Tortoise algorithm for directed sampling
#' uses a raster to find optimal sample locations.
#' The sampling strategy is based on the principle that an interpolation of the
#' samples should be as similar as possible to the guide raster.First, the center point
#' of the raster cell with the maximum deviation from the raster mean is sampled. Then,
#' the raster cell with the maximum deviation from the first sampled raster cell
#' is sampled. From then on, the values of the sampled raster cells are interpolated by
#' inverse distance weighting (idw), and the center point of the raster cell with
#' the largest absolute difference to the guide raster (the largest error) is sampled.
#' A new idw interpolation is made and a new cell is sampled. This is repeated
#' until a stopping criterion (stop_n or stop_dens) is reached.The sampling can be
#' stratified by a square grid. When a sample has been placed in a stratum, no more
#' samples will be placed in that stratum again until all other strata have been sampled.
#' The likelihood for a clipped stratum, at the edge of the area to be sampled,
#' is equal to the area of that stratum divided by the area of a full stratum.
#' Samples are placed in raster cell center points.
#'
#' The optional raster processing steps:
#' is carried out in the following order:
#' 1) mean filtering (argument: filter)
#' 2) resampling to specified resolution (argument: resolution),
#' 3) computation of first principal component (if x has multiple layers).

#' @import terra
#' @import gstat
#' @import sf
#' @importFrom 'stats' 'complete.cases'

#' @references Olsson, D. 2002. A method to optimize soil sampling from
#' ancillary data. Poster presenterad at: NJF seminar no. 336,
#' Implementation of Precision Farming in Practical Agriculture, 10-12 June
#' 2002, Skara, Sweden.

#' @examples
#' #load packages
#' require(terra)
#' require(gstat)
#' require(sf)
#' #create an example raster dataset
#' x<-rast(nrow=5, ncol=10, vals= sample(1:4, 50, replace=TRUE),crs=crs("EPSG:3857"))
#' x<-disagg(x, 10, 'bilinear')
#'
#' #create an example SpatVector of polygons
#' a<-cbind(
#'   x=c(-100, -120, -75, 40, 100, 120,  50, -100),
#'   y=c(-50, 0, 50, 75, 40, -30, -60, -50)
#' )
#' y<-vect(a, "polygons"); crs(y)<-crs(x)
#'
#' #do a directed stratified sampling for 25 samples. Let the stratification
#' #grid size be determined automatically and visualize the sampling
#' #procedure (default)
#' tortoise(x, y, stop_n=25)

#' @export
tortoise<-function(x=NULL, y=NULL, method='directed', edge=0, strat_size=NULL,
                   min_dist=0, p_idw=2, nmax_idw=8, resolution= NULL, filter=1,
                   stop_n=NULL, stop_dens=1, plot_results=T){

#check data ####
  #prepare for collecting feedback
  feedback<-'Messages'

  #are data provided?
  test<-method=='directed'& is.null(x)
  msg<-'x is missing'
  if(test) {stop(call. =F, msg); feedback<-c(feedback,msg)}

  test<-method=='directed'& is.null(y)
  msg<-'y is missing'
  if(test) {stop(call. =F, msg); feedback<-c(feedback,msg)}

  #are data projected?
  test<-method=='directed' & crs(x)==""
  msg<-'x is not projected.'
  if(test){stop(call. =F, msg); feedback<-c(feedback,msg)}

  test<-method=='directed' & crs(y)==""
  msg<-'y is not projected.'
  if(test){stop(call. =F, msg); feedback<-c(feedback,msg)}

  #are the data in a projected coordinate system?
  test<-is.lonlat(x)
  msg<-'The coordinate reference system of x is not projected.'
  if(test){stop(call. =F, msg); feedback<-c(feedback,msg)}

  test<-is.lonlat(y)
  msg<-'The coordinate reference system of y is not projected.'
  if(test){stop(call. =F, msg); feedback<-c(feedback,msg)}

  #are the raster and polygon data projected onto the same coordinate system?
  test<-method=='directed' & crs(x)!=crs(y)
  msg<-'x and y is have different coordinate reference systems.'
  if(test){stop(call.=F, msg); feedback<-c(feedback,msg)}

  #are the rasters of class numeric?
  test<-method =='directed' & !is.numeric(x[])
  msg<-'x is not numeric.'
  if(test){stop(call.=F, msg); feedback<-c(feedback,msg)}

  #is the polygon overlapping the raster?
  test<-method=='directed' & !is.related(x, y, 'intersects')
  msg<-'x any y do not intersect'
  if(test){stop(call. =F, msg); feedback<-c(feedback,msg)}

  #set stop_n to NULL if it is zero
  if(!is.null(stop_n))if(stop_n==0)stop_n<-NULL

  #is neither strat_size  nor stop_n has specified
  test1<-is.null(strat_size)
  test2<-T; if(!test1) test2<-strat_size==0
  test3<-is.null(stop_n)
  test4<-T; if(!test3) test4<-stop_n==0
  msg<-paste0('At least one of the arguments stop_n and strat_size needs ',
              'to be specified and larger than 0.')
  if((test1|test2)&(test3|test4)){stop(call. =F, msg); feedback<-c(feedback,msg)}

  #is resolution 0
  test1<-!is.null(resolution)
  test2<-resolution==0
  msg<-paste0('Resolution is 0. It must be NULL or larger than 0.')
  if(test1)if(test2){stop(call. =F, msg); feedback<-c(feedback,msg)}

  #is the buffer too wide
  test1<-edge>0
  test2<-sum(expanse(buffer(y,-edge)))==0
  msg<-paste0('The buffer is too wide, please chose a smaller value ',
              'of argument edge.')
  if(test1)if(test2){stop(call. =F, msg); feedback<-c(feedback,msg)}

  #is a too large number of samples per strata (stop_dens) specified?
  test<-stop_dens>10
  msg<-paste0('A large number of samples per strata has been specified. ',
              'The sampling may take a very long time. Consider changing argument ',
              'stop_dens to a smaller value or use the default.'
  )
  if(test){warning(call. =F, msg); feedback<-c(feedback,msg)}

  #is stop_dens larger than 1 and method='grid'
  test<-method=='grid'& stop_dens>1
  msg<-'Argument stop_dens was set to 1, which is the maximum allowed for method: grid.'
  if(test){stop_dens<-1; warning(call. =F, msg); feedback<-c(feedback,msg)}

  #is a too large total number of samples (stop_n) specified
  test1<-!is.null(stop_n)
  test2<-stop_n>500
  msg<-paste0('A large number of samples has been specified.',
              'The sampling may take a very long time. Consider ',
              'changing argument stop_n to a smaller value or use the default.')
  if(test1)if(test2){warning(call. =F, msg); feedback<-c(feedback,msg)}

  #warn if specified resolution is too fine
  a<-sum(expanse(y))
  finest_resolution<- signif((sqrt(a))/300, 1)
  test1<-method =='directed' & !is.null(resolution)
  test2<-resolution < finest_resolution
  msg<-paste0('Message: The specified resolution was very fine in relation to the size of ',
              'the area. Processing will take a very long time. Consider setting a coarser ',
              'resolution to speed up the processing.')
  if(test1)if(test2){
    message(appendLF = FALSE, msg); feedback<-c(feedback,msg)
  }

  #warn if specified resolution is too coarse
  test1<-method=='directed'
  test2<-ncell(x)<100
  msg<-paste0('Message: The raster resolution is relatively coarse in relation to the size of ',
              'the area. Consider resampling the raster to another resolution (argument: resolution).')
  if(test1) if(test2) {message(appendLF = FALSE, msg); feedback<-c(feedback,msg)}

  #shall stratification grid with an approximate strat_size fulfilling stop_dens and stop_n be created
  autostratification<-F
  test<-is.null(strat_size)
  msg<-paste0('Message: Argument strat_size is NULL The size of the ',
              'stratification grid is computed from the the total number of ',
              'samples (argument stop_n) and the numbers of samples per stratum ',
              '(argument stop_dens).')
  if(test){
    autostratification<-T
    message(appendLF = FALSE, msg); feedback<-c(feedback,msg)
  }

  test1<-!is.null(strat_size)
  test2<-strat_size==0 & method=='grid'
  msg<-paste0('Message: Argument strat_size is 0, whih is not valid for method grid.
              The size of the stratification grid is computed from the the total ',
              'number of samples (argument stop_n) and the numbers of samples per ',
              'stratum (argument stop_dens).')
  if(test1)if(test2){
    autostratification<-T
    message(appendLF = FALSE, msg); feedback<-c(feedback,msg)
  }

  test1<-method=='random' |method=='directed'
  test2<-!is.null(strat_size)
  test3<-strat_size==0
  if(test1)if(test2)if(test3) {autostratification<-F}

  test1<-!is.null(strat_size)
  test2<-!is.null(stop_n)
  test3<-strat_size!=0
  msg<-paste0('Message: Both strat_size and stop_n are specified. The size of the ',
              'stratification grid is computed from the the total number of ',
              'samples (argument stop_n) and the numbers of samples per stratum ',
              '(argument stop_dens).')
  if(test1 & test2)if(test3){
    autostratification<-T
    message(appendLF = FALSE, msg); feedback<-c(feedback,msg)
  }

  #is a too fine strat_size specified?
  test1<-!is.null(strat_size)
  test2<-strat_size!=0
  test3<-!autostratification & sum(expanse(y))/(strat_size^2)>500
  msg<-paste0('A very fine stratification grid is specified in relation to ',
              'the size of the area. The sampling will take very long time. Consider ',
              'changing argument strat_size to a larger value or use the default.')
  if(test1)if(test2)if(test3){warning(call. =F, msg); feedback<-c(feedback,msg)}

  #is a too fine stratification grid used compared to raster cell size?
  test1<-!is.null(strat_size)
  test2<-strat_size!=0
  test3<-method =='directed' & strat_size<res(x)[1]
  msg<-'The stratification grid has a smaller cell size than the raster.'
  if(test1)if(test2)if(test3){stop(call. =F, msg); feedback<-c(feedback,msg)}

  #is a too large filter provided
  test<-method=='directed' & filter>0.1*ncell(x)
  msg<-'Argument filter is too large in relation to the number of cells in x'
  if(test){stop(call. =F, msg); feedback<-c(feedback,msg)}

#prepare data, part1 ####
  #create polygon(y), if not provided
  test<-!is.null(x) & is.null(y)
  if(test){
    bbx<-as.polygons(ext(x))
    crs(bbx)<-crs(x)
    y<-bbx
  }

  #clip the polygon, if the raster does not fill it
  test<-method =='directed'
  if(test){
    bbx<-as.polygons(ext(x))
    crs(bbx)<-crs(x)
    y<-intersect(y, bbx)
  }

  #shrink polygon, if a buffer (argument edge) is specified
  ylarge<-y #save a copy
  y <- buffer(y,-edge)
  y<-y[expanse(y)>0,]  #repair y (remove polygons without area)

  #resample to the resolution specified by the user
  test1<-method =='directed' & !is.null(resolution)
  test2<-res(x)[1]!=resolution
  msg<-paste0('Message: Your raster was resampled to the specified resolution.')
  if (test1)if(test2){
    template<-rast(res=resolution, ext=ext(x), crs=crs(x))
    x<-resample (x=x, y=template, method='bilinear')
    message(appendLF = FALSE, msg); feedback<-c(feedback,msg)
  }

  #roughly crop and mask raster by the polygon
  test<-!is.null(x) & method!='grid' & method!='stratrand'
  if(test)x<-crop(x, buffer(ylarge, 5*res(x)[1]))

  #compute pca of rasterstack
  test1<-method=='directed'
  test2<-dim(x)[3]>1
  msg<- paste0('Message: More than one raster layer was provided so the first pricipal ',
               'component of them will be used to direct the sampling.')
  if(test1) if(test2){
    df <- as.data.frame(x, na.rm=T)
    sel<-complete.cases(x[])
    x$pc1<-NA
    x$pc1[sel]<-data.frame(stats::prcomp(df, scale = TRUE)['x'])[,1]
    x<-x$pc1
    message(appendLF = FALSE, msg); feedback<-c(feedback,msg)
  }

  #mean filter raster data
  test1<-method=='directed' & filter>1
  test2<-filter %% 2 == 0
  msg<-'Message: Your raster have been mean filtered. See argument filter.'
  if(test1){
    if(test2) filter<-filter+1 #if filter is even, add one
    x<-focal(x=x ,w=filter, fun='mean', na.policy= 'omit',na.rm=TRUE)
    message(appendLF = FALSE, msg); feedback<-c(feedback,msg)
  }

  #resample if raster resolution is too fine and no new resolution is specified
  test1<-method =='directed' & is.null(resolution)
  test2<-res(x)[1] < finest_resolution
  msg<-paste0('Message: Your raster was resampled to to speed up the processing. If you do not want this to happen, please provide the desired resolution (using argument resolution).')
  if(test1)if(test2){
    template<-rast(res=finest_resolution, ext=ext(x), crs=crs(x))
    x<-resample (x=x, y=template, method='bilinear')
    message(appendLF = FALSE, msg); feedback<-c(feedback,msg)
  }

  #final cropping and masking of raster(s) by polygon
  test<-!is.null(x) & method=='directed'
  if(test){
    x<-crop(x, ylarge)
    x<-mask(x=x, mask=y, touches=F)
  }

  #if needed, compute an approximate strat_size
  test<-autostratification
  if(test)strat_size<-signif(sqrt(sum(expanse(y)))/sqrt(stop_n/stop_dens), 2)
####

#create a stratification grid ####
  #define function for stratification
  f<-function(y, strat_size, stop_dens){
    strat<-rast(ext = ext(y), resolution = strat_size, crs= crs(y))
    strat<-extend(strat, c(1,1))
    strat<-as.polygons(strat)
    values(strat)<-1:ncell(strat)
    names(strat)<-'id'
    strat<-intersect(y, strat)
    strat$n<-0
    strat$target_n<-round(stop_dens*(expanse(strat)/max(expanse(strat))))
    strat$round_remain_n<-0
    strat$id<-1:nrow(strat)

    return(strat)
  }
  #create a stratification grid
  test1<-!is.null(strat_size)
  test2<-strat_size!=0
  if(test1)if(test2)strat<-f(y=y, strat_size=strat_size , stop_dens=stop_dens)

  #iteratively adjust strat_size to reach stop_n, if needed
  test1<-autostratification
  test2<-!is.null(stop_n)
  test3<-stop_n>0
  if(test1)if(test2)if(test3){
    for (i in 1:1000){
      strat<-f(y=y, strat_size=strat_size, stop_dens=1)
      n<-sum(strat$target_n)
      if(stop_n==n)break
      if(stop_n>n)strat_size<-strat_size*0.995
      if(stop_n<n)strat_size<-strat_size/0.995
    }
    #test if the iterative adjustment of strat_size was successful
    test<-stop_n!=n
    msg<-paste0('Could not find optimal stratificaiton size by iteration. ',
                'the number is samples will not be as specified by agrument stop_n. ',
                'Please try to specify an approximate stratification size (argument: ',
                'strat_size) to find a solution.')
    if(test){warning(call. =F, msg); feedback<-c(feedback,msg)}
  }

  #if a non-stratified method is used, just use the (buffered) polygon (one stratum)
  test1<-method=='random' |method=='directed'
  test2<-!is.null(strat_size)
  test3<-strat_size==0
  if(test1)if(test2)if(test3) {
    strat<-y
    strat$id<-1
    strat$n<-0
    strat$target_n<-stop_n
    strat$round_remain_n<-0
  }
####

#prepare data, part2 ####
  #set stop_dens to 1 if method = 'grid'
  if(method=='grid')stop_dens<-1

#identify number of samples to take, if stop_n is missing
  test<-is.null(stop_n)
  msg<-paste0('Message: The number of samples to place has been determined automatically.')
  if (test){
    stop_n<-sum(strat$target_n)
    message(appendLF = FALSE, msg); feedback<-c(feedback,msg)
  }
####

#sample ####
if(method=='directed'){
  #create a raster with cell numbers to use for masking
  cellnumbers<-x
  cellnumbers[]<-1:ncell(x)

  #rename raster
  names(x)<-'r'
  #run sampling loop
  for (i in 1:stop_n){

    #computes samples left to take in strata in total
    strat$tot_remain_n <-strat$target_n-strat$n

    #refill no of samples left to take in this round, if needed
    if(sum(strat$round_remain_n)==0) strat[strat$tot_remain_n>0, 'round_remain_n']<-1

    #calculate prediction raster (pred)
    if(i==1) {
      pred<- x
      pred[!is.na(pred)]<-mean(pred[], na.rm=T)
    } else{
      mod <- gstat(id = 'r', formula = r~1, locations=~x+y,
                   data=as.data.frame(p.sp, geom='XY'),
                   nmax=nmax_idw, set=list(idp=p_idw))
      pred <- interpolate(x, mod, debug.level = 0)[[1]]
      pred <- mask(pred, x)
    }

    #calculate aboslute error raster (ae)
    ae<-abs(x-pred)

    #save mae in p.sp
    if(i>1){p.sp$mae[i-1]<-global(ae, 'mean', na.rm=T)}

    #mask sampling area 1: keep only strata with samples left to take
    m<-strat[strat$round_remain_n==1,]

    #mask sampling area 2: exclude areas within min_distance to points
    if (nrow(m)>0 & i>1) {
      w<-max(min_dist, 0.5*res(x)[1])
      b<-buffer (p.sp, w)
      m<-m-b
      m<-m[expanse(m)>0,]
    }

    #test if there is any area left to sample
    test<-nrow(m)==0
    msg<-paste0('The specified number of samples could not be placed. No area ',
                  'left to sample. Try a smaller min_dist, a smaller edge, a ',
                  'smaller resolution, and/or a smaller ',
                  'stop_n!')
      if(test){stop(call. =F, msg)}

    #mask ae raster
    aemasked<-mask(ae, m, touches=F)
    if(i>1 & min_dist>0){
      min_dist_buffer<-y-buffer(x=p.sp, width=min_dist)
      aemasked<-mask(aemasked, min_dist_buffer, touches=F)
    }

    #test if there is any area left to sample
    test<-global(!is.na(ae), 'sum')==0
    msg<-paste0('The specified number of samples could not be placed. No area ',
                'left to sample. Try a smaller min_dist, a smaller edge, a ',
                'smaller resolution, and/or a smaller ',
                'stop_n!')
    if(test){stop(call. =F, msg)}

    #identify cell with maximum absolute error (ae), or for the first sample (i=1) the minimum ae
    if(i==1) {cell<-where.min(aemasked, values=F, list=T)[[1]]}
    if(i>1) {cell<-where.max(aemasked, values=F, list=T)[[1]]}
    cell<-cell[sample(x=length(cell), size=1)]

    test<-length(cell)==0
    msg<-paste0('The specified number of samples could not be placed. No area ',
                'left to sample. Try a smaller min_dist, a smaller edge, a ',
                'smaller resolution, and/or a smaller ',
                'stop_n!')
    if(test) {stop(call. =F, msg); feedback<-c(feedback,msg)}

    msk<-cellnumbers==cell
    msk<-subst(msk, FALSE, NA)
    pi.sp<-as.points(mask(x, msk))
    pi.sp$n<-i
    pi.sp$mae<-NA  #create column
    if(i==1)p.sp<-pi.sp
    if(i>1)p.sp<-rbind(p.sp, pi.sp)

    #update strat
    id<-strat[pi.sp]$id #identify current stratum
    sel<-strat$id==id
    strat$n[sel]<-strat$n[sel]+1
    strat$round_remain_n[sel]<-strat$round_remain_n[sel]-1

    #plot results
    if(plot_results==T){
      graphics::par(mfrow=c(2,2))
      plot(ylarge, main= 'Original surface')
      plot(x, legend=F, add=T,  )
      plot(strat, add=T)

      plot(ylarge, main= 'Reconstructed surface')
      plot(pred, legend=F, add=T)
      plot(strat, add=T)
      plot(p.sp, add=T, pch=20)

      plot(ylarge, main= 'Difference between the surfaces')
      plot((pred-x), add=T,legend=F)
      plot(strat, add=T)
      plot(p.sp, add=T, pch=20)

      if(i>3) {
        ymx<-max(pretty(p.sp$mae))
        plot(p.sp$n[3:nrow(p.sp)], p.sp$mae[3:nrow(p.sp)],
             ylab='Difference', xlab='No of samples',
             xlim=c(0, stop_n),  ylim=c(0, ymx),
             main= 'Difference vs no of samples'
        )
      }
    }
  }

}
if(method=='grid'){
  #find centroids of strata that are large enough to be sampled
  sel<-strat$target_n>0
  p.sp<-centroids(strat[sel,]) ##create centroids to be used for grid sampling

  #rename columns
  p.sp$n<-1:nrow(p.sp)
  p.sp<-p.sp[,'n']

  #convert to data.frame
  p.sp.df<-as.data.frame(p.sp, geom='XY')

  #plot results
  if(plot_results){
    graphics::par(mfrow=c(1,1))
    plot(ylarge, main='Grid sampling')
    plot(strat[strat$target_n>0], add=T, col='grey90')
    plot(p.sp, pch=20, add=T)
  }


}
if(method=='random'){#run sampling loop
  for (i in 1:stop_n){

    #computes samples left to take in strata in total
    strat$tot_remain_n <-strat$target_n-strat$n

    #refill no of samples left to take in this round, if needed
    if(sum (strat$round_remain_n)==0) strat[strat$tot_remain_n>0, 'round_remain_n']<-1

    #prepare mask for sampling area 1: keep only strata with samples left to take
    m<-strat[strat$round_remain_n==1,]

    #pepare mask sampling area 2: exclude areas within min_dist to points
    if(nrow(m)>0 & i>1 & min_dist > 0) {
      cl<-buffer (p.sp, min_dist)
      m <- m-cl
    }

    #test if there is any area left to sample
    test<-dim(m)[1]==0
    msg<-paste0('The specified number of samples could not be placed. No area ',
                'left to sample. Try a smaller min_dist, a smaller edge, a ',
                'smaller resolution, and/or a smaller ',
                'stop_n!')
    if(test) {stop(call. =F, msg); feedback<-c(feedback,msg)}

    #place a random sample
    m.sf<-st_as_sf(m)
    pi.sp<-st_sample(x=m.sf, size=1, type='random')
    pi.sp<-vect(pi.sp)
    pi.sp$n<-i
    if(i==1)p.sp<-pi.sp else p.sp<-rbind(p.sp, pi.sp)

    #update strat
    id<-strat[pi.sp,]$id
    sel<-strat$id==id
    strat[sel,]$n<-strat[sel,]$n+1
    strat[sel,]$round_remain_n<-strat[sel,]$round_remain_n-1

    #plot results
    if(plot_results==T){
      graphics::par(mfrow=c(1,1))
      plot(ylarge, main='Random sampling')
      plot(strat[strat$target_n>0], add=T, col='grey90')
      plot(p.sp, pch=20, add=T)
    }
  }

  #convert to data.frame
  p.sp.df<-as.data.frame(p.sp, geom='XY')
}
####

#compile results####
  #prepare  data
  sampled_area<-y
  stratification<-strat[,c('id', 'n')]
  samples<-p.sp[,'n']

  #list all results
  results<-list(sampled_area=sampled_area,
                  stratification=stratification,
                  samples=samples,
                  feedback=feedback)
  #add output if method is directed
  if (method=='directed'){
    sampled_raster<-x
    samples<-p.sp[,c('n', 'mae')]
    l<-list(sampled_raster=sampled_raster)
    results<-c(l, results)
  }
####

#return results
return(results)
}


