#' @title Calculate turnover, nestedness or total dissimilarity
#'
#' @param mpaShp shp file of the existing or proposed MPA
#' @param siteName character string of site name (no spaces) for file name creation (should match aqmClip value)
#' @param buff distance around mpa shape file for which analysis is conducted (to calculate extent box, "ploArea"). In decimal degrees (numeric)
#' @param thresh value from 0-1 to indicate presence threshold to filter data (usually one of 0.4, 0.5, 0.6, 0.7,default is 0.4)
#' @param dissCalc one of "turnover","nestedness","dissimilarity" indicating which metric of dissimilarity to plot. All are calculated based on the dissIndex each run
#' @param dissIndx one of "sorensen" or "jaccard" to indicate which metric of turnover and dissimilarity to use
#' @param corr indicates whether longitude "lonFx" needs to be corrected for crossing the 180 meridian, if no use NA
#'
#' @return ggplot map object with 0.5 degree raster of dissimilarity, with MPA shapefile, plot area, eez, and surrounding countries plotted
#' @export
#'
#' @examples
#' aqmTurnDiss(mpaShp=revMPA,siteName="rev",buff=10,thresh=0.6,dissCalc="turnover",dissIndex="jaccard",coor=NA)
aqmTurnDiss<-function(mpaShp,siteName,buff,thresh=0.4,dissCalc=c("turnover","nestedness","dissimilarity"),dissIndx=c("sorensen","jaccard"),corr=c(NA, "lonFx")){

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

  ### DATA
  # aqSpp<-read_csv("data/speciesoccursum.csv")

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

  ### format for biodiversity functions

  sppLst<-bind_rows(siteNow,siteFut)%>%
    select(SpeciesID)%>%
    filter(!duplicated(SpeciesID))%>%
    st_drop_geometry()

  turnDatN<-sppLst%>%
    left_join(.,siteNow)%>%
    # filter(Prbblty>=0.6)%>%
    mutate(coords=as.character(geometry),
           pres=case_when(probability >= thresh ~ 1,
                          TRUE ~ 0),
           pres=as.numeric(pres))%>%
    select(-probability)%>%
    pivot_wider(names_from = SpeciesID,values_from = pres)#%>%
  # select(-c(coords))

  turnDatN[is.na(turnDatN)] <- 0

  turnDatF<-sppLst%>%
    left_join(.,siteFut)%>%
    # filter(Prbblty>=0.6)%>%
    mutate(coords=as.character(geometry),
           pres=case_when(probability >= thresh ~ 1,
                          TRUE ~ 0),
           pres=as.numeric(pres))%>%
    select(-c(geometry,probability))%>%
    pivot_wider(names_from = SpeciesID,values_from = pres)#%>%
  # select(-c(coords))

  turnDatF[is.na(turnDatF)] <- 0

  ### calculate dissimilarity using the betapart package
  turnovrOut<-beta.temp(x=turnDatN%>%select(-c(geometry,coords)),
                        y=turnDatF%>%select(-c(coords)),
                        index.family=dissIndx)

  dissDat<-bind_cols(turnDatN%>%select(geometry),turnovrOut)%>%
    st_as_sf(crs = st_crs(mpa))%>%
    mutate(lon=st_coordinates(.)[,1],
           lat=st_coordinates(.)[,2])

  ### plots:
  ## get countries map

  ### crop base layers from world map

  if(is.na(corr)) {
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

  ## get eez

  # mrEezName<-mr_geo_code(place = mrSearch, like = TRUE, fuzzy = FALSE)%>%filter(placeType=="EEZ")%>%select(accepted)%>%as.numeric()
  #
  # siteEez<-mr_shp(key="MarineRegions:eez",filter=mrEezName)%>%
  #   st_as_sf(.,
  #            crs="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
  #
  # if(corr=="lonFx") {
  #   siteEez<-siteEez%>%
  #     st_as_sf(.,
  #              coords=c("long","lat"),
  #              crs=4326)%>%#"+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")%>%
  #     st_cast(to="POINT")%>% # convert polygons to points
  #     mutate(lon=st_coordinates(geometry)[,"X"] ,
  #            lat=st_coordinates(geometry)[,"Y"])%>%
  #     mutate(lon2=case_when(lon < 0 ~ 180+(180+lon),
  #                           TRUE ~ lon))%>%
  #     st_drop_geometry()%>%
  #     sfheaders::sf_multipolygon(x="lon2",y="lat",multipolygon_id=colnames(siteEez)[1])
  # }

  ## Plot layers

  datInpt<-if(dissIndx=="sorensen"){
    if(dissCalc=="turnover"){beta.sim} else {
      if(dissCalc =="dissimilarity"){beta.sor} else {
        beta.sne
      }
    }
  } else {
    if(dissCalc=="turnover") {beta.jtu} else {
      if(dissCalc =="dissimilarity") {beta.jac} else
      {beta.sne}
    }

    betaTitle<-if(dissCalc=="turnover"){"Species Turnover"} else {
      if(dissCalc =="dissimilarity"){"Species Dissimilarity"} else {
        "Species Nestedness"
      }
    }
  }

  dissplot<-ggplot()+
    geom_tile(data=dissDat,aes(x=lon,y=lat,fill=!!as.symbol(datInpt)),width=0.5,height=0.5)+
    scale_fill_gradient(high = "#da4325",low = "#e7e1bc",limits=c(0,1),name="")+
    geom_sf(data=siteEez,fill=NA,color="navy")+
    geom_map(data=geoShp, map=geoShp,aes(x=long, y=lat,map_id=region))+
    geom_sf(data = mpa,fill=NA, color="grey87")+
    scale_x_continuous(expand=c(0,0))+
    scale_y_continuous(expand=c(0,0))+
    labs(title=betaTitle, subtitle = "2019 to 2050",x="longitude",y="latitude")+
    theme_bw()

  return(dissplot)

}
