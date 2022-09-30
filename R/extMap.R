#' @title Extent Map
#' @description takes a MPA shapefile and creates an extent map showing the area assessed based on the established plot area
#'
#' @param mpaShp shp file of the existing or proposed MPA
#' @param eezShp shp file of the relevant EEZ as produced in mpaShpBbBase function (list item 4)
#' @param buff distance around mpa shape file for which analysis is conducted (to calculate extent box, "ploArea"). In decimal degrees (numeric)
#' @param zoomLat latitudinal buffer for the entire extent map in decimal degrees (size depends on context extent map scale)
#' @param zoomLon longitudinal buffer for the entire extent map in decimal degrees (size depends on context extent map scale)
#' @param corr indicates whether longitude "lonFx" needs to be corrected for crossing the 180 meridian, if no use NA
#'
#' @return extent map of the plot area around the mpaShp based on buff
#' @export
#'
#' @examples
#' extMap(mpaShp=revMPA,eezShp=mxMap[4][[1]], buff=7,zoomLat=15,zoomLon=10,corr=NA)

extMap<-function(mpaShp,eezShp,buff,zoomLat,zoomLon,corr=c(NA,"lonFx")){

  mapItems<-mpaShpBbBase(mpaShp = mpaShp,
                         buff=buff,
                         corr=corr)

  mpa=mapItems[1][[1]]
  ploArea=mapItems[2][[1]]
  geoShp=mapItems[3][[1]]

  siteBB<-st_bbox(mpa)

  bbDf<-data.frame(lon=c(as.numeric(ploArea$xmin),
                         as.numeric(ploArea$xmax),
                         as.numeric(ploArea$xmax),
                         as.numeric(ploArea$xmin)),
                   lat=c(as.numeric(ploArea$ymin),
                         as.numeric(ploArea$ymin),
                         as.numeric(ploArea$ymax),
                         as.numeric(ploArea$ymax)))%>%
    sfheaders::sf_polygon()%>%
    st_set_crs(4326)



  if(is.na(corr)) {
    totArea<-st_bbox(c(xmin=as.numeric(siteBB$xmin-zoomLon),
                       xmax=as.numeric(siteBB$xmax+zoomLon),
                       ymin=as.numeric(siteBB$ymin-zoomLat),
                       ymax=as.numeric(siteBB$ymax+zoomLat)),crs = st_crs(mpa))
  } else {
    totArea<-st_bbox(c(xmin=ifelse(as.numeric(siteBB$xmin)<0,
                                   180+(180+as.numeric(siteBB$xmin)),
                                   as.numeric(siteBB$xmin))-zoomLon,
                       xmax=ifelse(as.numeric(siteBB$xmax)<0,
                                   180+(180+as.numeric(siteBB$xmax)),
                                   as.numeric(siteBB$xmax))+zoomLon,
                       ymin=as.numeric(siteBB$ymin)-zoomLat,
                       ymax=as.numeric(siteBB$ymax)+zoomLat),
                     crs = st_crs(mpa))
  }

  ### crop base layers from world map
  if(is.na(corr)){

    geoShp<-map_data("world")%>%
      filter(long>totArea$xmin,
             long<totArea$xmax,
             lat>totArea$ymin,
             lat<totArea$ymax)%>%
      group_by(subregion)%>%
      mutate(dup=n()>3)%>%
      filter(dup==TRUE)%>%
      ungroup()

  } else {

    geoShp<-map_data("world",wrap=c(0,360))%>%
      filter(long>totArea$xmin,
             long<totArea$xmax,
             lat>totArea$ymin,
             lat<totArea$ymax)%>%
      group_by(subregion)%>%
      mutate(dup=n()>3)%>%
      filter(dup==TRUE)%>%
      ungroup()

  }


  extMap<-ggplot()+
    geom_map(data=geoShp,map=geoShp,aes(map_id=region))+
    # geom_sf(data=siteEez,fill="grey33",color=NA,alpha=0.3)+
    geom_sf(data=mpa,fill="deepskyblue3",alpha=0.3,color="deepskyblue4",linetype=2)+
    geom_sf(data=bbDf,fill=NA,color="darkorange3")+
    annotate(geom="text", x=ploArea$xmin+5, y=ploArea$ymin-1, label="Area assessed",
             color="darkorange3")+
    scale_x_continuous(expand = c(0,0))+
    scale_y_continuous(expand = c(0,0))+
    coord_sf(crs = 4326,
             xlim = c(totArea$xmin, totArea$xmax),
             ylim = c(totArea$ymin, totArea$ymax))+
    theme_bw()+
    labs(title = "Extent map",x="longitude",y="latitude")

  return(extMap)

}
