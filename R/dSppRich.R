#' @title Changes to species richness
#' @description Takes the output data from aqmClip and calculates the absolute change in species richness for the indicates region
#'
#' @param mpaShp shp file of the existing or proposed MPA, creates different layers using mpaShpBbBase function
#' @param siteName character string of site name (no spaces) for file name creation (should match aqmClip value)
#' @param buff distance around mpa shape file for which analysis is conducted (to calculate extent box, "ploArea"). In decimal degrees (numeric)
#' @param thresh value from 0-1 to indicate presence threshold to filter data (usually one of 0.4, 0.5, 0.6, 0.7,default is 0.4)
#' @param corr indicates whether longitude "lonFx" needs to be corrected for crossing the 180 meridian, if no use NA
#'
#' @return map of 0.5 degree raster of changes +/- in species richness
#' @export
#'
#' @examples
#' dSppRich(revMPA,siteName="rev",buff=10,thresh=0.6,corr=NA)
dSppRich<-function(mpaShp,siteName,buff,thresh=0.4,corr=c(NA, "lonFx")){

  ## pull in & format maps data
  if(is.na(corr)) {

    mpa<-mpaShp%>%
      st_set_crs(4326)#"+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

  } else {
    mpa<-mpaShp%>%
      st_set_crs(4326)%>%
      st_cast(to="POINT")%>% # convert polygons to points
      mutate(lon=st_coordinates(geometry)[,"X"] ,
             lat=st_coordinates(geometry)[,"Y"])%>%
      mutate(lon2=case_when(lon < 0 ~ 180+(180+lon),
                            TRUE ~ lon))%>%
      st_drop_geometry()%>%
      sfheaders::sf_multipolygon(x="lon2",y="lat",multipolygon_id=colnames(mpaShp)[1])%>%
      st_set_crs(4326)
  }

  siteBB<-st_bbox(mpa)

  if(is.na(corr)) {
    ploArea<-st_bbox(c(xmin=as.numeric(siteBB$xmin-buff),
                       xmax=as.numeric(siteBB$xmax+buff),
                       ymin=as.numeric(siteBB$ymin-buff),
                       ymax=as.numeric(siteBB$ymax+buff)),
                     crs = st_crs(mpa))
  } else {
    ploArea<-st_bbox(c(xmin=as.numeric(ifelse(siteBB$xmin<0,
                                              180+(180+siteBB$xmin),
                                              siteBB$xmin)-buff),
                       xmax=as.numeric(ifelse(siteBB$xmax<0,
                                              180+(180+siteBB$xmax),
                                              siteBB$xmax)+buff),
                       ymin=as.numeric(siteBB$ymin-buff),
                       ymax=as.numeric(siteBB$ymax+buff)),
                     crs = st_crs(mpa))
  }

  aqSpp<-read_csv("data/speciesoccursum.csv")

  ### current data
  nowFile<-paste("data/",siteName,"Aqnow",".shp",sep="")
  if (file.exists(nowFile)) {

    siteNow<-read_sf(nowFile)

  } else {
    stop(paste("The file does not exist: ",nowFile,". Create this file with 'aqmapClip()' before re-running"))
  }

  colnames(siteNow)<-c("SpeciesID","probability","geometry")

  ### future data
  futFile<-paste("data/",siteName,"Aqfut",".shp",sep="")

  if (file.exists(futFile)) {

    siteFut<-read_sf(futFile)

  } else {
    stop(paste("The ",futFile, " file does not exist. Create this file with 'aqmapClip()' before re-running"))
  }

  colnames(siteFut)<-c("SpeciesID","probability","geometry")

  ### Calculate species richness by cell (# of species)

  nowSR<-siteNow%>%
    filter(probability >= thresh)%>%
    mutate(geomID = as.character(geometry))%>%
    group_by(geomID)%>%
    summarise(sppRich=n())

  futSR<-siteFut%>%
    filter(probability >= thresh)%>%
    mutate(geomID=as.character(geometry))%>%
    group_by(geomID)%>%
    summarise(sppRich50=n())

  dSppRich<-nowSR%>%
    st_join(.,futSR)%>%
    mutate(srDiff=sppRich50-sppRich,
           lon=st_coordinates(.)[,1],
           lat=st_coordinates(.)[,2])

  ### plot changes in species richness

  ## get countries map

  ### crop base layers from world map

  if(is.na(corr)) {
    geoShp<-map_data("world")%>%
      st_as_sf(.,
               coords=c("long","lat"),
               crs="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")%>%
      group_by(group) %>%
      summarise(geometry = st_combine(geometry)) %>%
      st_cast("POLYGON")%>%
      st_crop(.,ploArea)

  } else {

    geoShp<-map_data("world")%>%
      st_as_sf(.,
               coords=c("long","lat"),
               crs=4326)%>%#"+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")%>%
      st_cast(to="POINT")%>% # convert polygons to points
      mutate(lon=st_coordinates(geometry)[,"X"] ,
             lat=st_coordinates(geometry)[,"Y"])%>%
      mutate(lon2=case_when(lon < 0 ~ 180+(180+lon),
                            TRUE ~ lon))%>%
      filter(lon2>ploArea$xmin & lon2<ploArea$xmax,
             lat>ploArea$ymin & lat<ploArea$ymax)%>%
      st_drop_geometry()%>%
      sfheaders::sf_multipolygon(x="lon2",y="lat",multipolygon_id="group")%>%
      st_make_valid()
  }

  ## get eez

  siteEez<-mapItems[4][[1]]

  ## Plot layers

  dSRplot<-ggplot()+
    geom_tile(data=dSppRich,aes(x=lon,y=lat,fill=srDiff),width=0.5,height=0.5)+
    scale_fill_gradient2(low="navyblue",mid="white",high="firebrick",midpoint = 0,name="Change in\nnumber of species")+
    geom_sf(data=siteEez,fill=NA,color="navy")+
    geom_sf(data = geoShp,fill="grey",color="black")+
    geom_sf(data = mpa,fill=NA, color="deepskyblue4",linetype=2)+
    scale_x_continuous(expand=c(0,0))+
    scale_y_continuous(expand=c(0,0))+
    labs(title="Changes in species richness (2019 to 2050)", subtitle = "Data source: Aquamaps",x="longitude",y="latitude")+
    theme_bw()

  return(dSRplot)


}
