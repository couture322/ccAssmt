#' @title Plot marine heatwave data spatially and temporally
#'
#' @param mpaShp shp file of the existing or proposed MPA
#' @param mhwMetric one of "dur" (mean duration), "intns" (maximal intensity), "sum" (number of MHW days per year)
#' @param ssps one of 126 or 585 indicating the projection scenario to look at
#' @param buff distance around mpa shape file for which analysis is conducted (to calculate extent box, "ploArea"). In decimal degrees (numeric)
#' @param YEAR year in the future for analysis from 2020 - 2100, default is 2050
#' @param corr indicates whether longitude (lonFx) needs to be corrected for crossing the 180 meridian, if no use NA
#' @param datFolder file path to where mhw data are stored
#'
#' @return two plots: 1) raster map of indicated MHW metric within indicated region, 2) plot of MHW metric over time
#' @export
#'
#' @examples
#' mhwMpaYrs(mpaShp=revMPA,mhwMetric="dur",ssps=585,buff=10,year=2075,corr=NA,datFolder="../../../ccAssmts/data")
mhwMapYrs<-function(mpaShp,mhwMetric,ssps,buff,YEAR=2050,corr=c(NA, "lonFx"),datFolder){

  ### format mpa shape file and set plot area bounding box

  mapItems<-mpaShpBbBase(mpaShp = mpaShp,
                         buff=buff,
                         corr=corr)

  mpa=mapItems[1][[1]]
  ploArea=mapItems[2][[1]]
  geoDf=mapItems[3][[1]]
  siteEez<-mapItems[4][[1]]

  #### pull MHW data from Google drive

  drFile<-if(ssps==126) {
    if(mhwMetric=="dur") {
      "multi_model_mean_ssp126_YEARMEANDURATION.nc"
    } else {if(mhwMetric=="intns") {
      "multi_model_mean_ssp126_YEARMAXINTENSITY.nc"
    } else {
      "multi_model_mean_ssp126_YEARSUM.nc"
    }}
  } else {
    if(mhwMetric=="dur") {
      "multi_model_mean_ssp585_YEARMEANDURATION.nc"
    } else {if(mhwMetric=="intns") {
      "multi_model_mean_ssp585_YEARMAXINTENSITY.nc"
    } else {
      "multi_model_mean_ssp585_YEARSUM.nc"
    }}
  }

  # datFile<-"mhwDat.nc" #where to put it
  #
  # drive_download(drFile, #raster file
  #                path=datFile,
  #                overwrite = TRUE)

  ### assuming I'm working on the server, directly pull in data rather than download from googledrive


  mhwDat<-terra::rast(paste(datFolder,drFile,sep = "/"))%>%#"mhwDat.nc")%>%
    as.data.frame(.,xy=TRUE)%>% ## convert to df
    rename(lon=x,
           lat=y)%>%
    pivot_longer(cols = tos_1:tos_251,
                 names_to = "yearCode",
                 values_to = "mhwVals")%>%
    mutate(yrIndx=as.numeric(sub("tos_","",yearCode)),
           year=1849+yrIndx)%>%
    filter(lon>ploArea$xmin,
           lon<ploArea$xmax,
           lat>ploArea$ymin,
           lat<ploArea$ymax)


  # ### delete large file once the object has been created
  # if (file.exists(datFile)) {
  #   unlink(datFile)
  #   cat(paste("The file has been deleted: ",datFile))
  # }


  ### plot local MHW projection map 2100

  mhwMetricTxt<-if(mhwMetric=="dur"){
    "Mean duration of MHW"

  } else {
    if(mhwMetric=="intns") {
      "Maximum MHW intensity"
    } else {
      "Total MHW days"
    }
  }


  if(mhwMetric=="intns"){
    unitsTxt<-"Degrees C"
    sfgLims<-c(0,ceiling(max(mhwDat%>%
                               filter(year==YEAR)%>%
                               select(mhwVals))))
  } else {
    unitsTxt<-"Days"
    sfgLims<-c(0,365)
  }


  mhwPlot<-mhwDat%>%
    filter(year==YEAR)%>%
    ggplot()+
    geom_raster(aes(x=lon,y=lat,fill=mhwVals))+
    scale_fill_gradient(high = "#da4325",low = "#e7e1bc",limits=sfgLims,name=unitsTxt)+
    geom_map(data=geoDf, map=geoDf,aes(x=long, y=lat,map_id=region))+
    geom_sf(data=mpa,fill=NA,alpha=0.25,color="grey87",size=0.5)+
    geom_sf(data=siteEez,fill=NA,color="navy")+
    scale_x_continuous(expand = c(0,0))+
    scale_y_continuous(expand = c(0,0))+
    theme_bw()+
    labs(title = mhwMetricTxt,subtitle=paste("Model projections for",as.character(YEAR),sep=" "),y="latitude",x="longitude")

  ### local temporal plot


  mhwYrsPlt<-ggplot(mhwDat)+
    geom_rect(xmin=2015,xmax=2100,ymin=-Inf,ymax=Inf,color=NA,fill="azure2",alpha=0.4)+
    geom_point(aes(x=year,y=mhwVals,color=lat))+
    scale_color_viridis_c()+
    scale_x_continuous(expand = c(0,0))+
    labs(title="Change over time",subtitle = mhwMetricTxt,y=unitsTxt,x="Year")+
    theme_bw()

  return(list(mhwPlot,mhwYrsPlt))

}
