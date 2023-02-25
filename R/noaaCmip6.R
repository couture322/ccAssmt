#' @title Plot NOAA CMIP6 data
#' @description Pulls NOAA data from author's private drive folder and plots up variable (temperature (tos), salinity (sos), oxygen at 2 depths (o2)) projections for the assessment area indicated
#'
#' @param mpaShp shp file of the existing or proposed MPA
#' @param eezShp shp file of relevant EEZ, if none enter NA
#' @param buff distance around mpa shape file for which analysis is conducted (to calculate extent box, "ploArea"). In decimal degrees (numeric)
#' @param futYr year to capture variable projection. Data range from 2015-2100
#' @param var one of "tos" (SST), "sos" (SSSal), "o2s" ([O2] shallow), "o2d" ([O2] deep) variables to plot.SST and SSS data are by definition only calculated at one depth (0m) and oxygen concentration is assessed at shallow and deep
#' @param corr indicates whether longitude "lonFx" needs to be corrected for crossing the 180 meridian, if no use NA
#'
#' @return ggplot object of variable projections or, for oxygen, list of 2 ggplot objects of projected oxygen concentations at 2 depths
#' @export
#' @import tidyverse sf googledrive patchwork
#' @examples
#' noaaCMIP6(mpaShp=revMPA,buff=7,var="tos",corr=NA)
#'
noaaCMIP6<-function(mpaShp,
                     eezShp,
                     buff,
                     futYr,
                     var=c("tos","sos","o2"),
                     corr=c(NA,"lonFx")) {


  ### create map layers
  mapItems<-mpaShpBbBase(mpaShp = mpaShp,
                         eezShp = eezShp,
                         buff=buff,
                         corr=corr)

  mpa=mapItems[1][[1]]
  ploArea=mapItems[2][[1]]
  geoShp=mapItems[3][[1]]
  siteEez<-mapItems[4][[1]]

  lyrIndx<-c(60:72,(((futYr-2015)*12)+1):(((futYr-2015)*12)+12))

  ### pull in global data from the drive based on indicated variable (var)

  driveFile<-if(var=="tos") {
    "tos_Omon_ACCESS-CM2_ssp126_r1i1p1f1_gn_201501-210012.nc"
  } else {
    if(var=="sos"){
      "sos_Omon_ACCESS-CM2_ssp126_r1i1p1f1_gn_201501-210012.nc"
    } else {
      "o2_Oyr_ACCESS-ESM1-5_ssp126_r1i1p1f1_gn_2015-2100.nc"
    }
  }

  # drive_download(driveFile,
  #                path=datFile,
  #                overwrite = TRUE)

  datFile<-paste("../../../../sandbox-sparc/cmip6_data/noaaCortadv6",driveFile,sep="/")

  # ncDat<-nc_open(datFile)
  noaaDat<-terra::rast(datFile)
  crs(noaaDat)<-"+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

  if(is.na(corr)){

    ext(noaaDat)<-c(-180,180,-90,90)


  } else {
    ext(noaaDat)<-c(0,360,-90,90)
  }


  if(var %in% c("tos","sos")){

    noaa2020<-subset(noaaDat,61:72)
    n20<-mean(noaa2020)

    noaaFut<-subset(noaaDat,(((futYr-2015)*12)+1):(((futYr-2015)*12)+12))
    nFut<-mean(noaaFut)

    noaa2<-nFut-n20

  } else {

    o2indx<-if(var=="o2s") {lyrIndx + 3} else {lyrIndx + 26}
    o2YrInx<-o2indx[c(6,futYr-2014)]
    noaaSub<-subset(noaaDat,o2YrInx) ## this should give us 2 raster layers: current (2020) and future (futYr) to calculate the change in [o2]

    noaa2<-subset(noaaSub,2)-subset(noaaSub,1)

  }

  noaaDiff<-noaa2%>%
    as.data.frame(.,xy=TRUE)%>% ## convert to df
    filter(x>ploArea$xmin,
           x<ploArea$xmax,
           y>ploArea$ymin,
           y<ploArea$ymax)
  colnames(noaaDiff)<-c("lon","lat","diff")

  ### plot data
  if(var=="tos") {
    pltTitl<-"Change in Sea Surface Temperature"
    pltSubTitl<-"Degrees Celcius" } else {
      if(var=="sos") {
        pltTitl<-"Change in Sea Surface Salinity"
        pltSubTitl<-"Parts per trillion"
      } else {
        pltTitl<-"Change in Oxygen Concerntration"
        pltSubTitl<-"Moles per cubic meter"
      }
    }

  if(any(noaaDiff$diff > 0) & any(noaaDiff$diff < 0)) {
    sclFllGrd<-scale_fill_gradient2(low = "#0a6165",high = "#da4325",mid = "#e7e1bc",midpoint = 0,
                                    name=paste("change\n2020 -", futYr))
  } else {

    if(all(noaaDiff$diff>=0)) {
      sclFllGrd<-scale_fill_gradient(high = "#da4325",low = "#e7e1bc",
                                     name=paste("change\n2020 -", futYr))
    } else {
      sclFllGrd<-scale_fill_gradient(low = "#0a6165",high = "#e7e1bc",
                                     name=paste("change\n2020 -", futYr))
    }
  }

  sitePlt<-ggplot()+
    geom_raster(data=noaaDiff,aes(x=lon,y=lat,fill=diff),width=1,height=1)+
    geom_sf(data=mpa,color="grey21",fill=NA)+
    geom_map(data=geoShp,map=geoShp,aes(map_id=region))+
    geom_sf(data=siteEez,color="navy",fill=NA)+
    sclFllGrd+
    scale_x_continuous(expand = c(0,0))+
    scale_y_continuous(expand = c(0,0),limits = c(ploArea$ymin,ploArea$ymax+0.1))+
    theme(legend.direction = "vertical",legend.position = "right")+
    labs(title=pltTitl,subtitle = pltSubTitl)+
    theme_bw()

  return(sitePlt)

}
