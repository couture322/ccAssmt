#' @title Spatial formatting
#' @description Takes a MPA shapefile and projects it for the rest of the analyses, creates a bounding box for the analysis based on buff and creates a base map of land masses within the bounding box
#'
#' @param mpaShp shp file of the existing or proposed MPA
#' @param buff distance around mpa shape file for which analysis is conducted (to calculate extent box, "ploArea"). In decimal degrees (numeric)
#' @param corr indicates whether longitude "lonFx" needs to be corrected for crossing the 180 meridian, if no use NA
#'
#' @return list of 3 items: mpaShp projected to WGS84 and corrected to 360 if needed: corr="lonFx", bounding box for the plot area, base map of land within the plot area
#' @export
#'
#' @examples
mpaShpBbBase<-function(mpaShp,buff,corr=c(NA,"lonFx")){

  ### format MPA shape file

  if(is.na(st_crs(mpaShp))) {

    mpa<-mpaShp%>%
      st_set_crs(4326)#"+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
  } else {

    mpa<-mpaShp%>%
      st_transform(4326)
  }

  if(is.na(corr)) {

    mpa<-mpa

  } else {

    mpa<-mpa%>%
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


  ### crop base layers from world map
  if(is.na(corr)){

    geoShp<-map_data("world")%>%
      filter(long>ploArea$xmin,
             long<ploArea$xmax,
             lat>ploArea$ymin,
             lat<ploArea$ymax)%>%
      group_by(subregion)%>%
      mutate(dup=n()>3)%>%
      filter(dup==TRUE)%>%
      ungroup()

  } else {

    geoShp<-map_data("world",wrap=c(0,360))%>%
      filter(long>ploArea$xmin,
             long<ploArea$xmax,
             lat>ploArea$ymin,
             lat<ploArea$ymax)%>%
      group_by(subregion)%>%
      mutate(dup=n()>3)%>%
      filter(dup==TRUE)%>%
      ungroup()

  }

  return(list(mpa,ploArea,geoShp))
}
