#' Calculates the great circle distances for a series of lon, lat coordinates
#' Returns adjusted coordinates to regress out of genetic data to correct for IBD
#'
#' @param lonlat two-column dataframe: longitude, latitude
#' @returns a list: \ gcd.coords - first two principal coordinates of the great circle distance matrix (two-column dataframe); \ gcd.distances - the distance matrix (dist object);
#' @export
gcd.dist=function(lonlat) {
  XY=lonlat
  gcd.slc <- function(twolonlats) {
    lon1=as.numeric(twolonlats[1])*pi/180
    lat1=as.numeric(twolonlats[2])*pi/180
    lon2=as.numeric(twolonlats[3])*pi/180
    lat2=as.numeric(twolonlats[4])*pi/180
    Re <- 6371 # Earth mean radius [km]
    d <- acos(sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2) * cos(lon2-lon1)) * Re
    return(d) # Distance in km
  }
  D=matrix(nrow=nrow(XY),ncol=nrow(XY))
  for (i in 1:(nrow(XY)-1)) {
    for (j in (i+1):nrow(XY)){
      if(sum(XY[i,]==XY[j,])==2) {
        D[i,j]=D[j,i]=0
      } else {
        D[i,j]=D[j,i]=gcd.slc(c(XY[i,],XY[j,]))
      }
    }
  }
  diag(D)=0
  gcdXY=data.frame(scores(capscale(D~1),scaling=1,display="sites",choices=c(1,2)))
  names(gcdXY)=c("x","y")
  return(list(gcd.coords=gcdXY,gcd.distances=as.dist(D)))
}


#' Calculates proportion of variation explained (R2) by each predictor, and in total
#'
#' takes into account the proportion of variation attributable to each of the ordination axes that went into the analysis.
#'
#' @param gf output of \code{gradientForest()}
#' @param ordination ordination the scores of which were analysed by \code{gradientForest()}
#' @param rescale2goods whether to rescale the variance attributable to the ordination's principal components to only those components that are pedictable by \code{gradientForest()}
#' @param rescale2npcs rescale total variation to the amount explained by this many leading PCs
#' @returns a named vector: proportion of total variation explained by each predictor
#' @export
importance_RDAforest=function(gf,ordination,rescale2goods=FALSE,rescale2npcs=NULL){
  mds.wts=(ordination$CA$eig/sum(ordination$CA$eig))
  if(rescale2goods==TRUE){
    mds.wts=mds.wts[names(gf$result)]
    mds.wts=mds.wts/sum(mds.wts)
  } else {
    if(!is.null(rescale2npcs)) {
      mds.wts=mds.wts[1:rescale2npcs]
      mds.wts=mds.wts/sum(mds.wts)
    }
  }
  gw=gf$imp.rsq
  for(m in colnames(gw)){
    gw[,m]=gw[,m]*mds.wts[m]
  }
  r2s=apply(gw,1,sum)
  r2s=r2s[order(r2s,decreasing=TRUE)]
  return(r2s)
}


#' Plot turnover curves for a gradientForest model
#'
#' wrapper for plot.gradientForest's plot.type="Cumulative.Importance"
#'
#' @param gf output of \code{gradientForest()}
#' @param vars variables to plot
#' @param show.MDS whether to show turnover curves for individual \code{Y} columns (assumed to be MDSes of an ordination)
#' @param common.scale whether to scale cumulative curves to the scale of the first plotted variable.
#' @param ... additional options to \code{plot.gradientForest()}
#' @import gradientForest
#' @export
plot_gf_turnovers=function(gf,vars,show.MDS=TRUE,common.scale=TRUE,...) {
  plot(gf,imp.vars=vars,plot.type="C",show.species=show.MDS,common.scale=common.scale,
       line.ylab = 0.9, par.args = list(
         mgp = c(1.5,0.5, 0),
         mar = c(2.5, 1, 0.1, 0.5),
         omi = c(0,0.3, 0, 0)
       )
  )
}

#' Lat-Lon to UTM (map coordinates) converter
#'
#' @param latlon lat, lon of points to convert
#' @param epsg epsg code for the region being plotted
#' @return dataframe with UTM-projected coordinates
#' @import sf
#' @import raster
#' @export
latlon2UTM=function(latlon,epsg){
  xy=latlon
  coordinates(xy) = ~lon+lat
  proj4string(xy) = CRS("+init=epsg:4326")
  targetCRS <- CRS(paste0("+init=epsg:",epsg))
  xy.p <- spTransform(xy, CRSobj = targetCRS)
  xy.map <- as.data.frame(xy.p)
  names(xy.map)=c("lon","lat")
  return(xy.map)
}

#' Determines epsg code for a give lon, lat
#'
#' @param lon longitude
#' @param lat latitude
#' @return epsg code
#' @export
epsg.maker=function(lon,lat){
   zone <- floor((lon + 180) / 6) + 1
   hemisphere_prefix <- ifelse(lat >= 0, 326, 327)
  epsg_code <- paste0(hemisphere_prefix, zone)
  return(as.numeric(epsg_code))
}


#' Plot a map in proper coordinates
#'
#' Plots two series of points on the map, background and overlay. Intended to plot colored raster (background) with sampled spots (overlay)
#'
#' @param coords latitutde and longitude of "background points" (colored raster). Must have columns named "lon" and "lat".
#' @param cols colors of background points
#' @param size size of background points. Increase until individual points merge into continuous colors. Default setting (0.5) seems to work with default Rmarkdown plot sizes (width=7, height=5)
#' @param pch symbol of background points.
#' @param margin map margin around the plotted region
#' @param mapdata map object of \code{sf} class, such as produced by \code{ne_countries} of the package \code{rnaturalearth}  (may need a pre-transform with \code{st_as_sf()})
#' @param map.on.top whether to plot the map on top of the raster points (preferable for seascape data)
#' @param overlay.points latitude and longitude of overlay points (in a typical case these would be sampling locations)
#' @param cols.points color(s) of overlay points
#' @param size.points size of overlay points
#' @param pch.points symbol for overlay points
#' @return (invisible) ggplot2 object
#' @import sf
#' @import ggplot2
#' @import raster
#' @export
plot_nice_map=function(coords,cols="black",mapdata,map.on.top=FALSE,overlay.points=NULL,cols.points=NULL,margin=40000,pch=16,size=0.5,pch.points=16,size.points=1) {
#  coords=xy2;mapdata=mapdata;map.on.top=F;margin=40000;cols=pa$colors;overlay.points = latlon;pch.points=10;size=1;pch=16;size.points=1;cols.points=NULL
  if(!is.null(overlay.points) & is.null(cols.points)) {
  # finding contrasting color for the point depending on background
    md=dist(coords[1:2,1:2])
    bgcol=list();i=1
    for(i in 1:nrow(overlay.points)) {
      table(coords[,1]>overlay.points[i,1])
      spot=which(
        coords[,"lon"]>overlay.points[i,"lon"]-1.5*md &
        coords[,"lon"]<overlay.points[i,"lon"]+1.5*md &
        coords[,"lat"]>overlay.points[i,"lat"]-1.5*md &
        coords[,"lat"]<overlay.points[i,"lat"]+1.5*md
      )
      bgcol[[i]]=apply(col2rgb(cols[spot]),1,mean)
    }
    bgcol.tab=do.call(rbind,bgcol)
#    cols.points0="hotpink"
    cols.points0=bw_choose(bgcol.tab[,1],bgcol.tab[,2],bgcol.tab[,3])
  }
  if(length(cols.points)==1) {cols.points0= cols.points}
  epsg.code=epsg.maker(mean(range(coords$lon)),mean(range(coords$lat)))
  mapdata_utm=st_transform(mapdata, crs = st_crs(epsg.code))
  xy.mapdata=latlon2UTM(coords,epsg.code)
  if(!is.null(overlay.points)) { xy.points=latlon2UTM(overlay.points,epsg.code) }
  bbox = st_as_sfc(st_bbox(c(xmin = min(xy.mapdata$lon)-margin, xmax = max(xy.mapdata$lon)+margin, ymin = min(xy.mapdata$lat)-margin, ymax = max(xy.mapdata$lat)+margin), crs = st_crs(epsg.code)))
  map.gg= geom_sf(data = mapdata_utm)
  raster.gg=geom_point(data = xy.mapdata, aes(lon,lat),col=adjustcolor(cols, alpha.f = 1),size=size,pch=pch)
  if(is.null(overlay.points)){
    points.gg=NULL
  } else {
    if(length(cols.points)>1) {
      cols.points_rev = rev(cols.points)
      points.gg=geom_point(data = xy.points, aes(lon,lat,color=cols.points_rev), pch=pch.points, size = size.points)
    } else {
      points.gg=geom_point(data = xy.points, aes(lon,lat),color=cols.points0, pch=pch.points, size = size.points)
    }
  }
  if(map.on.top==TRUE) {
    gg=ggplot()+raster.gg+map.gg+points.gg+
      coord_sf(xlim = c(st_bbox(bbox)["xmin"], st_bbox(bbox)["xmax"]),
               ylim = c(st_bbox(bbox)["ymin"], st_bbox(bbox)["ymax"]), expand = FALSE) +
      theme_minimal()
  } else {
    gg=ggplot()+map.gg+raster.gg+points.gg+
      coord_sf(xlim = c(st_bbox(bbox)["xmin"], st_bbox(bbox)["xmax"]),
               ylim = c(st_bbox(bbox)["ymin"], st_bbox(bbox)["ymax"]), expand = FALSE) +
      theme_minimal()
  }
  plot(gg)
  invisible(gg)
}

#' Dummify predictors
#'
#' Turns a dataframe containing numerical and categorical predictors into fully numerical
#'
#' @param X matrix of predictor variables
#' @param lastlevel whether to include all levels of a factor, or omit the last one (for lm designs)
#' @return a dataframe with 0,1 values and column names `factor1_level1`,`factor1_level2` etc (for one original column `factor1` with multiple levels). Originally numerical predictors are left unchanged.
#' @export
dummify=function(X,lastlevel=TRUE){
  Xd=list()
  for (ci in 1:ncol(X)) {
    if(is.factor(X[,ci]) | is.integer(X[,ci]) | is.character(X[,ci])) {
      zzzzz=as.factor(X[,ci])
      dum=data.frame(model.matrix(~0+zzzzz))
      if(lastlevel==FALSE) { dum=dum[,-ncol(dum)] }
      colnames(dum)=paste(colnames(X)[ci],colnames(dum),sep="_")
      colnames(dum)=sub("_zzzzz","_",colnames(dum))
      colnames(dum)=sub("^zzzzz","",colnames(dum))
      colnames(dum)=sub("_model.+","",colnames(dum))
      Xd[[ci]]=dum

    }  else {
      Xd[[ci]]=data.frame(X[,ci])
      names(Xd[[ci]])=colnames(X)[ci]
    }
  }
  Xdd=data.frame(do.call(cbind,Xd))
  return(Xdd)
}

#' Sum up importances
#'
#' Sums up importances of original factors that were dummified using `dummify()`
#'
#' @param gf \code{gradientForest()} result object
#' @param metadata original (non-dummified) metadata
#' @param ordination ordination the scores of which were analysed by \code{gradientForest()} (using function \code{makeGF()})
#' @param importance.cutoff do not sum up importances below that level
#' @return a dataframe of importances of original factors in `metadata`.
#' @export
sum_up_importances=function(gf,ordination=NULL,metadata,importance.cutoff=0){
  if (!is.null(ordination)) {
    ii=data.frame(importance_RDAforest(gf,ordination))
  } else {
    ii=data.frame(gradientForest::importance(gf,type="Weighted"))
  }
  names(ii)="importance"
  ii$var=row.names(ii)
  is=c();ii0=ii[ii$importance>=importance.cutoff,]
  for (f in colnames(metadata)) {
    i2=ii0[grep(paste("^",f,"_",sep=""),ii0$var),]
    if (nrow(i2)==0) {i2=ii0[grep(paste("^",f,sep=""),ii0$var),] }
    is=append(is,sum(i2$importance))
  }
  ii2=data.frame(cbind(importance=is,var=colnames(metadata)))
  ii2$importance=as.numeric(ii2$importance)
  ii2$var=factor(ii2$var,levels=ii2$var[order(ii2$importance)])
  return(ii2)
}


#' Run Gradient Forest (ordination version)
#'
#' This function runs gradient forest analysis on a vegan ordination object. It returns a gradientForest object.
#'
#' @param ordination ordination object made by \code{vegan::capscale()} or \code{vegan::rda()}
#' @param X matrix of predictor variables
#' @param ntrees number of random forest trees to create during run
#' @param pcs2keep principal axes of the ordination to analyze. Default is NULL, which means keeping them all.
#' @param ... additional parameters for the gradientForest() function
#' @import vegan
#' @import gradientForest
#' @export
makeGF=function(ordination,X,ntrees=1500,pcs2keep=NULL,...) {
  keep=pcs2keep
  Y=data.frame(scores(ordination,scaling=1,choices=c(1:length(ordination$CA$eig)))$sites)
  # removing constrained axes
  if(length(grep("CAP",colnames(Y)))>0) {Y=Y[,-grep("CAP",colnames(Y))] }
  if(!(is.null(keep))) { Y=Y[,keep] }
  if(dim(Y)[1]!=dim(X)[1]) { stop("incompatible input matrices: number of rows in X should be equal that in Y")}
  nSites = dim(Y)[1]
  nSpecs = dim(Y)[2]
  lev = floor(log2(nSites * 0.368/2))
  gf = gradientForest(cbind(X, Y),
                      predictor.vars = colnames(X), response.vars = colnames(Y),
                      ntree = ntrees, transform = NULL, trace=T,
                      maxLevel = lev, corr.threshold = 0.25,...)
  return(gf)
}

#' Run Gradient Forest (simple)
#'
#' This function is a simple wrapper for \code{gradientForest()} function, uses straight-up response matrix \code{Y}.
#'
#' @param Y matrix of response variables
#' @param X matrix of predictor variables
#' @param ntrees number of random forest trees to create during run
#' @param ... additional parameters for the \code{gradientForest()} function
#' @import vegan
#' @import gradientForest
#' @export
makeGF_simple=function(Y,X,ntrees=1500,...) {
  if(dim(Y)[1]!=dim(X)[1]) { stop("incompatible input matrices: number of rows in X should be equal that in Y")}
  nSites = dim(Y)[1]
  nSpecs = dim(Y)[2]
  lev = floor(log2(nSites * 0.368/2))
  gf = gradientForest(cbind(X, Y),
                      predictor.vars = colnames(X), response.vars = colnames(Y),
                      ntree = ntrees, transform = NULL, trace=T,
                      maxLevel = lev, corr.threshold = 0.25,...)
  return(gf)
}


#' Variable selection (jackknife version)
#'
#' Performs variable selection based on \code{mtry} criterion: variables that are not important by themselves but are correlated with actually important ones diminish in importance at higher \code{mtry} setting (Strobl et al 2018). The function runs ordination jackknife and selects variables that do not show decrease in importance at higher \code{mtry} in \code{prop.positive.cutoff} or more replicates.
#'
#' @param Y matrix of response variables or a distance matrix
#' @param X matrix of predictor variables
#' @param mintry lower \code{mtry} value to use; if omitted, will use \code{0.2*N} (where \code{N} is the number of predictors) or 3, whatever is greater.
#' @param maxtry higher \code{mtry} value to use; if omitted, will use \code{0.4*N+1} (where \code{N} is the number of predictors) or 6, whatever is greater.
#' @param oob out-of-bag fraction: the fraction of datapoints to withhold from ordination construction in each replicate
#' @param covariates data frame of covariates to be regressed out during ordination
#' @param nreps number of jackknifing replicates
#' @param ntrees number of trees to construct in each replicate
#' @param prop.positive.cutoff proportion of replicates that must show non-negative importance change at higher \code{mtry} to keep the variable.
#' @param importance.cutoff do not accept predictors with R2 less than this fraction of the R2 of the best predictor. For example, if the best predictor has \code{R2 = 0.2} and \code{importance.cutoff} is set to 0.1, predictors with R2 less than 0.02 will be discarded.
#' @param top.pcs number of top principal components to consider
#' @param importance.type by default, the method uses "GF", which is \code{gradientForest}-style "Weighted" importances (average R2 across predictable PCs). "GF" equally rewards variables explaining any PC. Alternative is "RDAF", the \code{RDAforest}-style R2: proportion of total variance explained by the predictor, adjusting  for the amount of total variance attributable to each PC. I feel "GF" is more appropriate for variable selection.
#' @param rescale2npcs (applicable when \code{importance.type="RDAF"}) whether to rescale total variation to the amount explained by the top PCs.
#' @return \code{goodvars} - variables passing the selection criteria;   \code{delta} - change in importance at higher mtry for each variable in each replicate;   \code{importances} - median importance (proportion of total variance explained) of each variable at higher \code{mtry};   \code{importances.GFscale} - median importances scaled in gradientForest style (mean R2 across predictable PCs);   \code{prop.positive} - proportion of positive importance changes at higher \code{mtry} for each variable
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr %>%
#' @import vegan
#' @import gradientForest
#' @export
mtrySelJack=function(Y,X,mintry=NULL,maxtry=NULL,nreps=15,oob=0.2,prop.positive.cutoff=0.333,importance.cutoff=0.1,importance.type="GF",top.pcs=10,rescale2npcs=TRUE,covariates=NULL,ntrees=1500) {
  #  Y=sc.fit;X=env.fit;covariates=NULL;nreps=2;oob=0.1;prop.positive.cutoff=0.5;top.pcs=3;ntrees=1500
    if(dim(Y)[1]!=dim(X)[1]) { stop("incompatible input matrices: number of rows in X should be equal that in Y")}
  imps=c();impsgf=c();replicate=c();delta=list()
  if(rescale2npcs==TRUE){
    rescale2npcs=top.pcs
  } else {rescale2npcs=NULL }
  # computing total variance (global inertia) in the top PCs
  if(dim(Y)[1]==dim(Y)[2] & sum(Y[,1]==Y[1,])==dim(Y)[1]) {
    if(!is.null(covariates)) {
      ords.full=capscale(Y~1+Condition(as.matrix(covariates)))
    } else {
      ords.full=capscale(Y~1)
    }
  } else {
    if(!is.null(covariates)) {
      ords.full=rda(Y~1+Condition(as.matrix(covariates)))
    } else {
      ords.full=rda(Y~1)
    }
  }
  failed.reps=0;v.order=NULL
  for(i in 1:nreps){
    message("\nreplicate ",i)
    ib.keep=sample(1:nrow(Y),round(nrow(Y)*(1-oob)))
    if(!is.null(covariates)) {
      covars.ib=covariates[ib.keep,]
      if(dim(Y)[1]==dim(Y)[2] & sum(Y[,1]==Y[1,])==dim(Y)[1]) {
        Y.ib=Y[ib.keep,ib.keep]
        ords0=capscale(Y.ib~1+Condition(as.matrix(covars.ib)))
        suppressWarnings({ords=predict(ords0,Y,type='sp',scaling="sites")})
      } else {
        Y.ib=Y[ib.keep,]
        ords0=rda(Y.ib~1+Condition(as.matrix(covars.ib)))
#        suppressWarnings({ords=predict(ords0,Y,type='wa',scaling="sites")})
        ords=predict(ords0,Y,type='wa',scaling="sites")
      }
    } else {
      if(dim(Y)[1]==dim(Y)[2] & sum(Y[,1]==Y[1,])==dim(Y)[1]) {
        Y.ib=Y[ib.keep,ib.keep]
        ords0=capscale(Y.ib~1)
        suppressWarnings({ords=predict(ords0,Y,type='sp',scaling="sites")})
#        ords=predict(ords0,Y,type='wa',scaling="sites")
      } else {
        Y.ib=Y[ib.keep,]
        ords0=rda(Y.ib~1)
#        suppressWarnings({ords=predict(ords0,Y,type='wa',scaling="sites")})
        ords=predict(ords0,Y,type='wa',scaling="sites")
      }
    }
    if(is.null(mintry)) { mintry=max(3,round(ncol(X)/5)) }
    if(is.null(maxtry)) { maxtry=max(min(ncol(X)-2,6),round(2*ncol(X)/5)+1) }
    # mintry=min(3,1+round(ncol(X)/4))
    # maxtry=min(6,2+round(0.375*ncol(X)))
    message("mtry ",mintry,"...")
    gf00=makeGF_simple(ords[,c(1:top.pcs)],X,ntrees=ntrees,mtry=mintry)
    message("\nmtry ",maxtry,"...")
    gf1=makeGF_simple(ords[,c(1:top.pcs)],X,ntrees=ntrees,mtry=maxtry)
    if(is.null(gf1) | is.null(gf00)) {
      failed.reps=failed.reps+1
      next
    }
    if(importance.type=="RDAF"){
      ii1=importance_RDAforest(gf1,ords0,rescale2npcs=rescale2npcs)
      ii0=importance_RDAforest(gf00,ords0,rescale2npcs=rescale2npcs)[names(ii1)]
    } else {
      ii1=gradientForest::importance(gf1,type="Weighted")
      ii0=gradientForest::importance(gf00,type="Weighted")[names(ii1)]
    }
    if (is.null(v.order)) {
      v.order=names(ii1)
    }
    delta[[i]]=ii1[v.order]-ii0[v.order]
    imps=c(imps,ii1[v.order])
  }
  # differences in R2-scaled importances
  delta=data.frame(do.call(cbind,delta))
  sdd=stack(delta)
  sdd$var=v.order
  sdd$var=factor(sdd$var,levels=rev(v.order))

  dimps=data.frame(imps)
  dimps$variable=names(imps)
  dimps$variable=factor(dimps$variable,levels=rev(v.order))
  meds=dimps%>%
    group_by(variable)%>%
    summarise(median(imps))
  meds=as.data.frame(meds)
  importances=meds$`median(imps)`
  names(importances)=meds$variable
  importances=importances[order(importances,decreasing=TRUE)]

  # proportion of positive importance changes
  dc=data.frame(apply(delta,1,function(x){return(sum(x>=0)/length(x))}))
  dc$var=v.order
  dc$var=factor(dc$var,levels=rev(v.order))
  names(dc)[1]="prop.positive"
  ic=importance.cutoff*max(importances)
  goodimps=names(importances)[importances>=ic]
  gooddelta=row.names(dc)[dc$prop.positive>=prop.positive.cutoff]
  goodvars=intersect(goodimps,gooddelta)
  return(list(goodvars=goodvars,delta=sdd,importances=importances,prop.positive=dc,failed.reps=failed.reps))
}

#' Predict multi-column matrix by random forest
#'
#' builds \code{gradientForest::randomForest} model for each column in \code{Y} based on \code{X}
#' then predicts \code{Y} column values based on \code{newX} values of predictors
#'
#' @param Y data frame of columns to predict
#' @param X matrix of predictor variables to build random forest models
#' @param newX new matrix of predictor variable where \code{Y} columns must be predicted.
#' @param extra if newX variable is outside the range of its X values by not more than this fraction of the variable's span, it is set to max or min of X. Values further away are set to NA and these rows are deleted.
#' @param ... other parameters for \code{gradientForest::randomForest}
#' @return list of two items: preds - dataframe of predicted \code{Y} values; goodrows - vector of passing (in-range) rows (T/F, subsetter for plotting)
#' @importFrom extendedForest randomForest
#' @export
predict_rf=function(Y, X, newX, extra=0, ...) {
  #  Y=pcs;X=env[,evars];newX=rasters.inrange;extra=0
  # making sure there are no values in newX that are outside the range in X
  newX=newX[,colnames(newX) %in% colnames(X)]
  ranges=apply(X,2,range)
  spans=ranges[2,]-ranges[1,]
  ranges.s=ranges
  ranges.s[2,]=ranges[2,]+spans*extra
  ranges.s[1,]=ranges[1,]-spans*extra
  for (v in colnames(newX)){
    newX[,v][newX[,v]<ranges.s[1,v]]=NA
    newX[,v][newX[,v]<ranges[1,v]]=ranges[1,v]
    newX[,v][newX[,v]>ranges.s[2,v]]=NA
    newX[,v][newX[,v]>ranges[2,v]]=ranges[2,v]
  }
  goodrows=apply(newX,1,function(x){sum(is.na(x))==0})
  newX=newX[goodrows,]
  Yp=list()
  for(i in 1:ncol(Y)){
    message("predicting mds ",i)
    pc=Y[,i]
    rf=randomForest(pc~.,data=cbind(pc,X),...)
    Yp[[i]]=predict(rf,newX)
  }
  return(list(preds=data.frame(do.call(cbind,Yp)),goodrows=goodrows))
}

#' Predict multi-column matrix by gradient forest (make turnover curves for clustering)
#'
#' builds \code{gradientForest::randomForest} model for each column in \code{Y} based on \code{X}
#' then predicts \code{Y} column values based on \code{newX} values of predictors
#'
#' @param Y data frame of columns to predict
#' @param X matrix of predictor variables to build random forest models
#' @param newX new matrix of predictor variable where \code{Y} columns must be predicted.
#' @param extra if newX variable is outside the range of its X values by not more than this fraction of the variable's span, it is set to max or min of X. Values further away are set to NA and these rows are deleted.
#' @param ... other parameters for \code{gradientForest::randomForest}
#' @return list of two items: preds - dataframe of predicted \code{Y} values; goodrows - vector of passing (in-range) rows (T/F, subsetter for plotting)
#' @import gradientForest
#' @export
predict_gf=function(Y,X,newX,extra=0,...){
#  Y=pcs;X=env[,evars];newX=rasters
  # making sure there are no values in newX that are outside the range in X
  newX=newX[,colnames(newX) %in% colnames(X)]
  ranges=apply(X,2,range)
  spans=ranges[2,]-ranges[1,]
  ranges.s=ranges
  ranges.s[2,]=ranges[2,]+spans*extra
  ranges.s[1,]=ranges[1,]-spans*extra
  for (v in colnames(newX)){
    newX[,v][newX[,v]<ranges.s[1,v]]=NA
    newX[,v][newX[,v]<ranges[1,v]]=ranges[1,v]
    newX[,v][newX[,v]>ranges.s[2,v]]=NA
    newX[,v][newX[,v]>ranges[2,v]]=ranges[2,v]
  }
  goodrows=apply(newX,1,function(x){sum(is.na(x))==0})
  newX=newX[goodrows,]
  gf=makeGF_simple(Y,X)
  pp=predict(gf,newX)
  return(list(preds=pp,goodrows=goodrows))
}


#' Ordination Jackknife
#'
#' Runs gradient forest on \code{nreps} replicates, each time re-building the ordination based on a fraction of the data, projecting the rest of datapoints.
#'
#' @param Y matrix of response variables or a distance matrix
#' @param X matrix of predictor variables
#' @param oob out-of-bag fraction: the fraction of datapoints to withhold from ordination construction in each replicate
#' @param newX new matrix of predictor variable where associated differences in Y must be predicted (averaged across replicates). If omitted, predictions will be made for the original X.
#' @param covariates data frame of covariates to be regressed out during ordination
#' @param extra if newX variable is outside the range of its X values by not more than this fraction of the variable's span, it is set to max or min of X. Values further away are set to NA and these rows are deleted.
#' @param nreps number of spatial bootstrap replicates
#' @param top.pcs number of top principal components to consider
#' @param rescale2npcs whether to rescale total variation to the amount explained by the top PCs.
#' @param ntrees number of trees to construct in each replicate
#' @param mtry \code{mtry} setting. If left unspecified, defaults to \code{N/3} (where \code{N} is the number of predictors in X)
#' @return A list of four items: \code{predictions.direct} - averaged predicted Y values for \code{newX}, \code{predictions.turnover} - averaged turnover curves for \code{newX}, \code{median.importance} - median importance of each variable, \code{all.importances} - all observed importances across all replicates. Table with three columns: \code{importance}, \code{variable}, \code{rep}.
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr %>%
#' @import vegan
#' @import gradientForest
#' @export
ordinationJackknife=function(Y,X,oob=0.2,newX=NULL,covariates=NULL,nreps=25,top.pcs=15,rescale2npcs=FALSE,mtry=NULL,ntrees=1500,extra=0) {
#    ntrees=1500;Y=IBS.p;oob=0.2;X=env.p[,mm$goodvars];newX=cbind(XY,rasters.post);nreps=3;rescale2npcs=FALSE;extra=0.1;covariates=cbind(latlon.p.utm,admix.cov);top.pcs=10;mtry=NULL
  if(dim(Y)[1]!=dim(X)[1]) { stop("incompatible input matrices: number of rows in X should be equal that in Y")}
  if(rescale2npcs==TRUE){rescale2npcs=top.pcs } else {rescale2npcs=NULL }
  # setting mtry to N/3 or 2 (note: the default is N/3 or 1)
  if (is.null(mtry)){ mt=max(floor(ncol(X)/3),2) } else { mt=mtry }
    failed.ords0=NULL
    failed.ords=NULL
    failed.reps=0;v.order=NULL
  if (!is.null(newX)) {
    newX=newX[,colnames(newX) %in% colnames(X)]
    ranges=apply(X,2,range)
    spans=ranges[2,]-ranges[1,]
    ranges.s=ranges
    ranges.s[2,]=ranges[2,]+spans*extra
    ranges.s[1,]=ranges[1,]-spans*extra
    for (v in colnames(newX)){
      newX[,v][newX[,v]<ranges.s[1,v]]=NA
      newX[,v][newX[,v]<ranges[1,v]]=ranges[1,v]
      newX[,v][newX[,v]>ranges.s[2,v]]=NA
      newX[,v][newX[,v]>ranges[2,v]]=ranges[2,v]
    }
    goodrows=apply(newX,1,function(x){sum(is.na(x))==0})
    newX=newX[goodrows,]
  } else {
    newX=X
    goodrows=c(1:nrow(X))
    }
  imps=c();replicate=c();impsum=c()

  if(dim(Y)[1]==dim(Y)[2] & sum(Y[,1]==Y[1,])==dim(Y)[1]) {
    if(!is.null(covariates)) {
      ords.full=capscale(Y~1+Condition(as.matrix(covariates)))
    } else {
      ords.full=capscale(Y~1)
    }
  } else {
    if(!is.null(covariates)) {
      ords.full=rda(Y~1+Condition(as.matrix(covariates)))
    } else {
      ords.full=rda(Y~1)
    }
  }

  sc.full=scores(ords.full,scaling=1,display="sites",choices=c(1:length(ords.full$CA$eig)))
  for(i in 1:nreps){
    message("\nreplicate ",i)
    ib.keep=sample(1:nrow(Y),round(nrow(Y)*(1-oob)))
    if(!is.null(covariates)) {
      covars.ib=covariates[ib.keep,]
      if(dim(Y)[1]==dim(Y)[2] & sum(Y[,1]==Y[1,])==dim(Y)[1]) {
        Y.ib=Y[ib.keep,ib.keep]
        ords0=capscale(Y.ib~1+Condition(as.matrix(covars.ib)))
        suppressWarnings({ords=predict(ords0,Y,type='sp',scaling="sites")})
      } else {
        Y.ib=Y[ib.keep,]
        ords0=rda(Y.ib~1+Condition(as.matrix(covars.ib)))
        suppressWarnings({ords=predict(ords0,Y,type='wa',scaling="sites")})
      }
    } else {
      if(dim(Y)[1]==dim(Y)[2] & sum(Y[,1]==Y[1,])==dim(Y)[1]) {
        Y.ib=Y[ib.keep,ib.keep]
        ords0=capscale(Y.ib~1)
        suppressWarnings({ords=predict(ords0,Y,type='sp',scaling="sites")})
      } else {
        Y.ib=Y[ib.keep,]
        ords0=rda(Y.ib~1)
        suppressWarnings({ords=predict(ords0,Y,type='wa',scaling="sites")})
      }
    }
    # align projected ordinaton with full-model scores
    orr=data.frame(procrustes(sc.full,ords,scale = T)$Yrot)
    names(orr)[1:top.pcs]=colnames(sc.full)[1:top.pcs]
#    message("names orr: ", names(orr)[1:top.pcs])
    gf1=makeGF_simple(orr[,1:top.pcs],X,ntrees=ntrees,mtry=mt)
    if(is.null(gf1)) {
      failed.reps=failed.reps+1
      message("   ... nothing is predictable. Moving to next rep")
      failed.ords0=ords0
      failed.ords=ords
      next
    }
    else{
      if (is.null(v.order)) {
        pred.d=predict_rf(orr[,1:top.pcs],X,newX,ntree=ntrees,mtry=mt)$preds
        pred.t=predict(gf1,newX)
        v.order=names(importance_RDAforest(gf1,ords0))
  #      message("v.order: ",v.order)
      } else {
        pred.d=pred.d+predict_rf(orr[,1:top.pcs],X,newX,ntree=ntrees,mtry=mt)$preds
        pred.t=pred.t+predict(gf1,newX)
      }
      ii=importance_RDAforest(gf1,ords0,rescale2npcs=rescale2npcs)[v.order]
  #    print(ii)
      imps=c(imps,ii)
      replicate=c(replicate,rep(i,length(v.order)))
      impsum=c(impsum,sum(ii))
    }
  }
  pred.d=pred.d/(nreps-failed.reps)
  pred.t=pred.t/(nreps-failed.reps)

  # computing median importances of variables
  dimps=data.frame(imps)
  dimps$replicate=replicate
  dimps$variable=names(imps)
  names(dimps)[1]="importance"
  meds=dimps%>%
    group_by(variable)%>%
    summarise(median(importance))
  meds=as.data.frame(meds)
  importances=meds$`median(importance)`
  names(importances)=meds$variable

  # rescaling turnover curves to actual proportion of variation explained
  for(m in colnames(pred.t)){
    pred.t[,m]=decostand(pred.t[,m],method="range")*importances[m]
  }

  # releveling variables according to median importance
  varorder=meds$variable[order(meds[,2],decreasing=T)]
#  turnovers=turnovers[,varorder[varorder %in% colnames(turnovers)]]
  importances=importances[varorder]
  dimps$variable=factor(dimps$variable,levels=rev(varorder))
  return(list(sum.importances=impsum,predictions.direct=pred.d,predictions.turnover=pred.t,median.importance=importances,all.importances=dimps,goodrows=goodrows,failed.reps=failed.reps,failed.ords0=failed.ords0,failed.ords=failed.ords))
}

#' Plot turnover curves returned by \code{ordinationJackknife}
#'
#' @param oj object returned by \code{ordinationJackknife}
#' @param X matrix of predictor variables used to form predictions in \code{oj}
#' @param variable name of the variable to plot
#' @param ... additional options for \code{plot()}
#' @export
plot_turnover=function(oj,X,variable,...) {
#  oj=sbll2e;X=predictors[sbll2e$goodrows,];variable="lon"
  tt=data.frame(cbind(X[,variable],oj$predictions.turnover[,variable]))
  ttt=tt[order(tt[,1]),]
  plot(ttt[,2]~ttt[,1],type="step",xlab="X value",ylab="cumulative R2",main=variable,mgp=c(2.3,1,0),...)
}

#' Plot adaptive neighborhoods
#'
#' Plots a map of differential local adaptation, colored by the first two or three PCs. Can cluster points with similar predicted adaptation, using either the predictions themselves or an additional dataset that is particularly good for clustering (such as predictions from \code{gradientForest}). In the latter case, clusters that contain similar predicted adaptation can be merged.
#' Plots three plots:
#' - hierarchical clustering tree for adaptation clusters with red line for the \code{cluster.merge.level} setting (clusters below this line will be merged);
#' - scatterplot of first two adaptation PCs with \code{envfit} arrows based on the supplied \code{envs} table of environmental variables;
#' - physical map of adaptive neighborhoods, based on supplied spatial coordinates (\code{XY}).
#'
#' @param Y matrix the PCs of which will be represented by colors on map. Within RDAforest pipeline, it is typically direct predictions from \code{ordinationJackknife}: \code{modelName$predictions.direct}.
#' @param envs environmental variables to show as \code{envfit} arrows
#' @param xy spatial coordinates of samples in \code{Y}
#' @param nclust number of color clusters to make
#' @param cluster.guide dataset to derive clusters from. Within RDAforest pipeline, it is typically turnover curves predicted by \code{ordinationJackknife}: \code{modelName$predictions.turnover}.
#' @param cluster.merge.level merge clusters that are more similar than this fraction of max cluster distance. Examine the hierarchical clustering tree plotted by this function to set this threshold.
#' @param cluster.labels whether to add cluster numbers to map
#' @param col.envs color of arrows and their text labels
#' @param ordinate to plot PCs of supplied Y matrix instead of Y directly (beta! use at your own risk)
#' @param matching.scores reorder and flip Y columns to match these scores (beta! use at your own risk)
#' @param color.scheme string of 0s or 1s of length 2 (to display the first 2 PCs) or of length 3 (for 3 PCs), for example "01" or "110". "1" goes for PCs that are going to be color-flipped.
#' @param lighten increase to 0.3-0.5 for less saturated colors, decrease to 0 for max saturation
#' @param scal scaling of envfit arrows: the higher, the shorter the arrows
#' @param jitscale controls distance from arrow point to text label
#' @param rangeExp controls overall X,Y range of the PCA plot (increase to make PCA itself more compact)
#' @param cex symbol size control
#' @param cex.txt text size control
#' @param ... other \code{plot} options for the final plot (the map)
#' @return (invisible) list of six items: \code{xy} - spatial coordinates; \code{clusters} - clustering indices; \code{PCs} - PC1 and PC2 coordinates; \code{vars} - predictor variable names; \code{arrows} - coordinates of (scaled) arrows on PCA plot, \code{colors} - colors for points.
#' @import clv
#' @import vegan
#' @export
plot_adaptation=function(Y,envs,xy,cluster.guide=NULL,matching.scores=NULL,ordinate=TRUE,nclust=0,cluster.merge.level=0.333,color.scheme="000",scal=15,jitscale=0.03,col.envs="black",rangeExp=2,cex=1,cex.txt=1,cluster.labels=TRUE,lighten=0.25,...){
#Y=dirs;envs=ras2[,bests];xy=xy2;lighten=0.5;cluster.guide=NULL;ordinate=F;color.scheme="000";cluster.merge.level=0.5;nclust=10;flip2=-1;flip1=1;scal=8;cex=1;jitscale=0.05;rangeExp=1.5;col.envs="black";cluster.labels=TRUE;pcs2show=3
  flips=as.numeric(strsplit(color.scheme,split="")[[1]])
  if(ordinate==TRUE){
      if(dim(Y)[1]==dim(Y)[2] & sum(Y[,1]==Y[1,])==dim(Y)[1]){ ordination=capscale(Y~1)} else { ordination=rda(Y~1)}
    nc=min(c(length(ordination$CA$eig),length(flips)))
    rd.sc=data.frame(scores(ordination,display="sites",scaling=1,choices=c(1:length(ordination$CA$eig))))
  } else {
    nc=min(c(ncol(Y),length(flips)))
    rd.sc=Y
  }
  if(!is.null(matching.scores)){
#    matching.scores=cpcs0[sb$goodrows,];rd.sc=dirs
    dm=round(cor(cbind(rd.sc,matching.scores)))
    diag(dm)=0;mm=c();ms=c();ir=1
    for(ir in 1:ncol(rd.sc)){
        m0=which(abs(dm[ir,])==1)
        mm=c(mm,m0)
        ms=c(ms,dm[ir,m0])
    }
    mm=mm-ncol(rd.sc)
    rd.sc=rd.sc[,mm]
    ms=ms[mm]
    for(ri in 1:length(ms)) {
      rd.sc[,ri]=rd.sc[,ri]*ms[ri]
    }
  }
  if(!length(flips) %in% c(2,3)){
    message("color.scheme appears to be not 2 or 3 characters long, switching to default \"000\"")
    flips=c(0,0,0)
  }
  rd.sc=rd.sc[,1:nc]
  flips[flips>0]=-1
  flips[flips==0]=1
  flipsm=diag(flips[1:nc])
#  print(flips)
#  diag(flips)=c(flip1,flip2,flip3)[1:nc]
  envft=envfit(rd.sc,envs,choices=c(1:nc))
  env.arrows=data.frame(scores(envft,display="vectors",scaling=1,choices=c(1:nc)))
  bests=row.names(env.arrows)
  for(i in 1:nc){ rd.sc[,i]=rd.sc[,i]*flips[i] }
#  rd.sc=as.matrix(rd.sc) %*% flips
  env.arrows=as.matrix(env.arrows) %*% flipsm

  pc1 = rd.sc[, 1]
  pc2 = rd.sc[, 2]
  if(nc==2) {
    b <- pc2
    g <- pc1-pc2
    r <- pc1
    # b <- pc1-pc2
    # g <- -pc1
    # r <- pc2
  } else {
    pc3=rd.sc[, 3]
    r <- pc1
    g <- pc3
    b <- pc2
  # r <- pc1+pc2
    # g <- -pc3
    # b <- pc3+pc2-pc1
  }
  lighten=lighten*255
  r <- (r - min(r))/(max(r) - min(r)) * 255
  g <- (g - min(g))/(max(g) - min(g)) * 255
  b <- (b - min(b))/(max(b) - min(b)) * 255
  clusters=NULL
  if(nclust==0) {
    r <- (r - min(r)+lighten)/(max(r) - min(r)+lighten) * 255
    g <- (g - min(g)+lighten)/(max(g) - min(g)+lighten) * 255
    b <- (b - min(b)+lighten)/(max(b) - min(b)+lighten) * 255
    colors=rgb(r,g,b,max=255)
    meds=NULL
    clusters=NULL
  } else {
    if (nclust>0) {
      if(is.null(cluster.guide)) {
        if(ordinate==TRUE) {
          PCs=data.frame(scores(ordination,display="sites",scaling=1,choices=c(1:length(ordination$CA$eig))))
        } else {
          PCs=Y
        }
      } else {
        ordi2=rda(cluster.guide~1)
        PCs=data.frame(scores(ordi2,display="sites",scaling=1,choices=c(1:length(ordi2$CA$eig))))
      }
      # message("computing optimal number of clusters...")
      # set.seed(1234)
      # mcc=Mclust(PCs,G=2:nclust)
      # set.seed(NULL)
      # message(mcc$G)
      clPCs = clara(PCs, nclust, sampsize = 1000)
      clusters=clPCs$clustering
#      clusters=as.integer(mcc$classification)
#      meds=GDAtools::medoids(dist(PCs),clusters)
      meds=clPCs$i.med
      cdd=cls.scatt.data(Y, clusters)$intercls.centroid
      colnames(cdd)=row.names(cdd)=names(meds)=unique(clusters)
      hc=hclust(as.dist(cdd),method="ave")
      mhi=hc$height[which(hc$height==max(hc$height))]*cluster.merge.level
      if(cluster.merge.level>0) {
        cl.indices.merged=cutree(hc,h=mhi)
        names(cl.indices.merged)=unique(clusters)
        clusters2=clusters
        clist=sort(unique(clusters))
        names(cl.indices.merged)=clist
        newnames=clist
         seen=c();i=1
        for(i in 1:length(clist)){
           cl=clist[i]
            clusters2[clusters==cl]=cl.indices.merged[cl]
             if(cl.indices.merged[cl] %in% seen) {
               newnames[i]=names(seen)[seen %in% cl.indices.merged[cl]]
             } else {
                 seen=c(seen,cl.indices.merged[cl])
               }
        }
#        names(meds)=cl.indices.merged
        names(meds)=newnames
        clusters=clusters2
        colnames(cdd)=row.names(cdd)=newnames
      }
      hc=hclust(as.dist(cdd),method="ave")
      hc$height=hc$height/hc$height[which(hc$height==max(hc$height))]
      plot(hc)
      abline(h=cluster.merge.level,col="red")
    }
    medcolR = c()
    medcolG = c()
    medcolB = c()
    clis=clusters
    for(i in 1:length(unique(clusters))){
      ci = unique(clusters)[i]
      medcolR=c(medcolR,mean(r[clusters==ci]))
      medcolG=c(medcolG,mean(g[clusters==ci]))
      medcolB=c(medcolB,mean(b[clusters==ci]))
      clis[clusters==ci]=i
    }
    medcolR <- (medcolR - min(medcolR)+lighten)/(max(medcolR) - min(medcolR)+lighten) * 255
    medcolG <- (medcolG - min(medcolG)+lighten)/(max(medcolG) - min(medcolG)+lighten) * 255
    medcolB <- (medcolB - min(medcolB)+lighten)/(max(medcolB) - min(medcolB)+lighten) * 255
    colors=rgb(medcolR[clis], medcolG[clis], medcolB[clis],max=255)
    bwc=bw_choose(medcolR[clis], medcolG[clis], medcolB[clis])
  }
  xrng <- range(rd.sc[, 1], env.arrows[, 1]/scal) * rangeExp
  yrng <- range(rd.sc[, 2], env.arrows[, 2]/scal) * rangeExp
  jitscale=jitscale*(xrng[2]-xrng[1])
  plot((rd.sc[, 1:2]), xlim = xrng, ylim = yrng, cex = cex, col = colors, pch=16,asp = 1,axes=FALSE,xlab="",ylab="")
  arrows(rep(0, length(bests)), rep(0,length(bests)), env.arrows[bests,1]/scal, env.arrows[bests,2]/scal, length = 0.0625,col=col.envs)
  jit = rnorm(length(bests),jitscale,jitscale/3)
  text(env.arrows[bests,1]/scal + jit * sign(env.arrows[bests,1]), env.arrows[bests,2]/scal + jit * sign(env.arrows[bests,2]),col=col.envs, labels = bests,cex=cex.txt)
  plot(xy, pch=16,cex = cex, asp = 1, col = colors,...)
  if(length(clusters)>1 & cluster.labels==TRUE) {
    text(xy[meds,], labels = names(meds),col=bwc[meds],cex=cex)
  }
  invisible(list(xy=xy,clusters=clusters,PCs=rd.sc[,1:2],vars=bests,arrows=env.arrows[bests,1]/scal,colors=colors,medoids=meds))
}


#' Black or white text labels, depending on background color?
#'
#' Calculates background luminance and returns black hex value for lighter backgrounds, white hex value for darker backgrounds.
#'
#' @param r vector of background red values
#' @param g vector of background green values
#' @param b vector of background blue value
#' @return vector of black or white hex values
#' @export
bw_choose=function(r,g,b){
  r <- (r - min(r))/(max(r) - min(r))
  g <- (g - min(g))/(max(g) - min(g))
  b <- (b - min(b))/(max(b) - min(b))
  RGBs=data.frame(cbind(r,g,b))
  bw= apply(RGBs,1,
            function(x){
              Y=0.2126*x[1] + 0.7152*x[2] + 0.0722*x[3]
              if(Y> 0.5) { return("#000000")} else { return("#FFFFFF")}
            }
  )
  return(bw)
}

#' Genetic offset
#'
#' Returns row by row difference between two dataframes containing
#' random forest predictions of adaptation
#'
#' @param X first dataframe
#' @param Y second dataframe (must be identical in shape to the first one)
#' @param method distance measure, see \code{dist}
#' @return vector of row-by-row distances
#' @export
gen_offset=function(X,Y,method="euclidean"){
  if((nrow(X) == nrow(Y)) & (ncol(X) == ncol(Y))) {
    di=c()
    for (i in 1:nrow(X)) {
      dd=data.frame(rbind(X[i,],Y[i,]))
      di=c(di,dist(dd,method=method))
    }
    return(di)
  } else {
    stop("gen_offset error: rasters are not the same shape")
  }
}

#' Genetic offset for ordinationJackknife objects
#'
#' Aligns objects, executes \code{gen_offset}, rescales to specified scale
#'
#' @param X present-day ordinationJackknife object
#' @param Y future ordinationJackknife object
#' @param sx coordinates of prediction space for X
#' @param sy coordinates of prediction space for Y
#' @param method distance measure, see \code{dist}
#' @param sc scaling factor for the offset. Options besides "none" are "q50" and "q90":  0.5th and 0.9th quantile of the disparity of present-day predictions across landscape.
#' @return data frame of coordinates and offsets
#' @export
gen_offset_oj=function(X,Y,sx,sy,sc="q90",method="euclidean"){
# X=oj;Y=ojf;sx=envc[,1:2];sy=envf[,1:2];sc="q90";method="euclidean"
  # aligning predictions
 names(sx)=names(sy)
  xp=data.frame(cbind(sx[X$goodrows,],X$predictions.direct))
  yp=data.frame(cbind(sy[Y$goodrows,],Y$predictions.direct))
  xyp=merge(xp,yp,by=names(sx))
  xym=xyp[,1:2]
  xm=xyp[,c(3:(ncol(X$predictions.direct)+2))]
  ym=xyp[,c((ncol(X$predictions.direct)+3):ncol(xyp))]
  names(xm)=names(ym)=names(X$predictions.direct)
# determining scale
  if (sc == "none") { sc=1 }
  if (sc == "q90") {
    sc=adapt_scale(X$predictions.direct)[2]
  } else {
      if (sc == "q50") {
        sc=adapt_scale(X$predictions.direct)[1]
      }
  }
  offs=gen_offset(xm,ym,method=method)/sc
  return(data.frame(cbind(xym,offset=offs)))
}

#' Environmental mismatch
#'
#' Assesses how similar a given gPC vector X is to gPC predictions across landscape.
#'
#' @param X vector of gPCs to compare to landscape (for example, predicted for a focal location, or for a specific individual)
#' @param Y \code{ordinationJackknife()} object containing gPCs predictions across landscape
#' @param sy coordinates of prediction space (for Y)
#' @param sc scaling factor. Suggested values are 0.5-th or 0.9-th quantile of present-day adaptation disparity across landscape, see \code{adapt_scale}
#' @param method distance measure, see \code{dist}
#' @return data frame of coordinates and distances from X
#' @export
env_mismatch=function(X,Y,sy,sc=1,method="euclidean"){
  dd=c()
  for(i in 1:nrow(Y$predictions.direct)){
    dd=c(dd,dist(rbind(Y$predictions.direct[i,],X),method=method))
  }
  dd=dd/sc
  return(data.frame(cbind(sy[Y$goodrows,],env.mismatch=dd)))
}


#' Adaptation scale
#'
#' Estimates quantiles of absolute Eucludean distances of random forest predictions
#' (to rescale genetic offsets)
#'
#' @param X dataframe of RF predictions
#' @param quantiles vector of quantiles to return
#' @param method distance measure, see \code{dist}
#' @export
adapt_scale=function(X,quantiles=c(0.5,0.9),method="euclidean"){
#X=oj$predictions.direct;method="euclidean";quantiles=c(0.5,0.9)
  nsamp=min(100,nrow(X))
  nrolls=3*round(nrow(X)/nsamp)
  dd=c()
  for (i in 1:nrolls){
    dd=c(dd,dist(X[sample(c(1:nrow(X)),nsamp),],method=method))
  }
  qq=quantile(abs(dd),quantiles,na.rm=TRUE)
  return(qq)
}


#' Reselect variables
#'
#' resets goodvars slot in mtrySelJack object based on new criteria
#'
#' @param M mtrySelJack object
#' @param prop.positive.cutoff cutoff for proportion of replicates where importance increases at higher mtry
#' @param importance.cutoff do not accept predictors with R2 less than this fraction of the R2 of the best predictor.
#' @return mtrySelJack object with modified `goodvars` slot
#' @export
Reselect=function(M,prop.positive.cutoff=0.5,importance.cutoff=0.1){
  vars=row.names(M$prop.positive)
  M$goodvars=vars[M$prop.positive[vars,"prop.positive"]>=prop.positive.cutoff & M$importances[vars]>importance.cutoff*max(M$importances)]
  return(M)
}



#' Create a multilayer SpatRaster out of lonlat and values dataframes
#'
#' ..to extract data with coordinates from it: as.data.frame(RASTERname,xy=TRUE)
#' 
#' 
#' @param lonlat dataframe of latitutde and longitude corresponding to rows in `values`. Must have columns named "lon" and "lat".
#' @param values dataframe of environmental variables
#' @param res desired resolution, in degrees latitude
#' @param polygon dataframe of points (lon, lat) defining a polygon. The first row must be the same as the last row.
#' @param variogram_model fixed variogram model, made with variogram() and fit.variogram(). If omitted, will be fit on the fly.
#' @param grid.res grid resolution, in meters. 
#' @param smooth.range in meters
#' @param utm.zone UTM zone number
#' @param verbose whether to plot original data, smoothed data, and kriging result
#' @return dataframe of lat-lon (grid within polygon), and projected variable
#' @import terra
#' @export
lonlat2raster=function(lonlat,values,res=0.01){
  names(lonlat)=c("lon","lat")
  dff=data.frame(cbind(lonlat,values))
  lon_range <- range(lonlat$lon, na.rm = TRUE)
  lat_range <- range(lonlat$lat, na.rm = TRUE)
  
  # Create an empty raster template
  template_raster <- rast(
    xmin = lon_range[1], xmax = lon_range[2],
    ymin = lat_range[1], ymax = lat_range[2],
    resolution = res,
    crs = "EPSG:4326"  # WGS84
  )
  
  # Convert dataframe to SpatVector
  data_vect <- vect(dff, geom = c("lon", "lat"), crs = "EPSG:4326")
  
  # Rasterize each layer
  layers <- lapply(names(values), function(layer_name) {
    rasterize(data_vect, template_raster, field = layer_name)
  })
  
  # Combine all layers into a multi-layer raster
  mrast <- rast(layers)
  
  # Assign names to layers
  names(mrast) <- names(values)
  
  return(mrast)
}


viridis_colors <- function(input_vector,...) {
  normalized_values <- scales::rescale(input_vector,to=c(0,1))  # Normalize to [0,1]
  colors <- viridis(100,...)  # Generate 100 Viridis colors for smooth mapping
  color_indices <- round(normalized_values * (length(colors) - 1)) + 1  # Scale to color indices
  return(colors[color_indices])  # Assign colors based on values
}



