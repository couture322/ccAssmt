#' @title Spatial formatting
#' @description Takes a MPA shapefile and projects it for the rest of the analyses, creates a bounding box for the analysis based on buff and creates a base map of land masses within the bounding box
#'
#' @param mpaShp shp file of the existing or proposed MPA
#' @param eezShp shp file of relevant EEZ to be converted to crs 4326, if NA, the EEZ will be pulled from the marine regions dataset *note: if the mpa or plot area cross the 180 meridian, marine regions data cannot be used
#' @param buff distance around mpa shape file for which analysis is conducted (to calculate extent box, "ploArea"). In decimal degrees (numeric)
#' @param corr indicates whether longitude "lonFx" needs to be corrected for crossing the 180 meridian, if no use NA
#'
#' @return list of 3 items: mpaShp projected to WGS84 and corrected to 360 if needed: corr="lonFx", bounding box for the plot area, base map of land within the plot area, eez projected to WGS84 and corrected to plot with the other elements
#' @export
#' @import tidyverse sf
#' @examples
#' mapsEx<-mpaShpBbBase(mpaShp=revMPA,eezShp=NA,buff=10,corr=NA)
#'
#' exMap<-ggplot()+geom_map(data=mapsEx[[3]],map = mapsEx[[3]],aes(map_id=region))+
#' geom_sf(data=mapsEx[[1]])+
#' geom_sf(data=mapsEx[[4]])
#'
#' exMap
#'
mpaShpBbBase<-function(mpaShp,eezShp,buff,corr=NA){

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
      dplyr::ungroup()

  } else {

    geoShp<-map_data("world",wrap=c(0,360))%>%
      filter(long>ploArea$xmin,
             long<ploArea$xmax,
             lat>ploArea$ymin,
             lat<ploArea$ymax)%>%
      group_by(subregion)%>%
      mutate(dup=n()>3)%>%
      filter(dup==TRUE)%>%
      dplyr::ungroup()

  }

  if(is.na(eezShp)){ ### an NA in the eezShp argument implies that corr = NA as well

    if(is.na(corr)){

      eez<-mrEEZ%>% ### mrEEZ data are already in WGS84
        st_crop(ploArea)

    } else {
      print("Marine regions data cannot be Longitude corrected, provide a shapefile for the EEZ if the assessment crosses the 180 meridian")
    }

  } else {

    if(is.na(st_crs(eezShp))) {

      eez<-eezShp%>%
        st_set_crs(4326)#"+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
    } else {

      eez<-eezShp%>%
        st_transform(4326)
    }

    if(is.na(corr)) {

      eez<-eez

    } else {

      eez<-eez%>%
        st_cast(to="POINT")%>% # convert polygons to points
        mutate(lon=st_coordinates(geometry)[,"X"] ,
               lat=st_coordinates(geometry)[,"Y"])%>%
        mutate(lon2=case_when(lon < 0 ~ 180+(180+lon),
                              TRUE ~ lon))%>%
        st_drop_geometry()%>%
        sfheaders::sf_multipolygon(x="lon2",y="lat",multipolygon_id=colnames(eezShp)[1])%>%
        st_set_crs(4326)
    }

  }

  return(list(mpa,ploArea,geoShp,eez))
}
