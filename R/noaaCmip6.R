
#' @title Plot NOAA CMIP6 data
#' @description Pulls NOAA data from author's private drive folder and plots up variable (temperature (tos), salinity (sos), oxygen at 2 depths (o2)) projections for the assessment area indicated
#'
#' @param mpaShp shp file of the existing or proposed MPA
#' @param buff distance around mpa shape file for which analysis is conducted (to calculate extent box, "ploArea"). In decimal degrees (numeric)
#' @param var one of "tos" (SST), "sos" (SSSal), "o2" ([oxygen]) variables to plot.SST and SSS data are by definition only calculated at one depth (0m) and oxygen concentration is assessed at shallow and deep
#' @param corr indicates whether longitude "lonFx" needs to be corrected for crossing the 180 meridian, if no use NA
#'
#' @return ggplot object of variable projections or, for oxygen, list of 2 ggplot objects of projected oxygen concentations at 2 depths
#' @export
#'
#' @examples
noaaCMIP6<-function(mpaShp,
                    buff,
                    var=c("tos","sos","o2"),
                    corr=c(NA,"lonFx")) {

  ### pull in global data from the drive based on indicated variable (var)
  datFile<-"data/tempDat/noaaTmp.nc"

  driveFile<-if(var=="tos") {
    "tos_Omon_ACCESS-CM2_ssp126_r1i1p1f1_gn_201501-210012.nc"
  } else {
    if(var=="sos"){
      "sos_Omon_ACCESS-CM2_ssp126_r1i1p1f1_gn_201501-210012.nc"
    } else {
      "o2_Oyr_ACCESS-ESM1-5_ssp126_r1i1p1f1_gn_2015-2100.nc"
    }
  }

  drive_download(driveFile,
                 path=datFile,
                 overwrite = TRUE)

  # ncDat<-nc_open(datFile)
  rastDat<-terra::rast(datFile)


  ### pull variables from nc file: ncDat
  lon<-ncvar_get(ncDat,"longitude")
  lat<-ncvar_get(ncDat,"latitude")
  datVar<-ncvar_get(ncDat,var) # pull variable data
  timeBnd<-ncvar_get(ncDat,"time_bnds")

  ttl<-ncatt_get(ncDat,var,"long_name")
  ttlNm<-as.character(ttl$value)

  ### delete large file from local

  if (file.exists(datFile)) {
    unlink(datFile)
    cat(paste("The file has been deleted: ",datFile))
  }

  ### convert arrays to DF

  if(var %in% c("tos","sos")){

    newYr<-seq(1,ncol(timeBnd),by=12)
    aveDf<-NULL
    datesDf<-tibble(date=as_date(timeBnd[1,],origin="1850-01-01"))%>%
      mutate(year=year(date))


    for(i in 1:length(newYr)){

      yrDf<-expand.grid(lat=lat[1,],lon=lon[,1])%>%
        mutate(lon=ifelse(lon>180,-(360-lon),lon))


      for(j in 1:12){ # iterate over time
        k=newYr[i]+j-1

        monDf <- expand.grid(lat=lat[1,],lon=lon[,1])%>%
          mutate(lon=ifelse(lon>180,-(360-lon),lon),
                 varVal=as.vector(t(datVar[,,k])))

        yrDf = yrDf %>% bind_cols(monDf%>%select(varVal))

      }
      yrDf<-yrDf%>%
        mutate(aveVar=rowMeans(.[-1:-2],na.rm = TRUE))%>%
        select(lat,lon,aveVar)%>%
        mutate(year=2014+i)#datesDf[k,"year"])

      aveDf<-aveDf%>%
        bind_rows(.,yrDf)

      i=i+1
    }
  } else {

    o2Dfshallow<-NULL

    datesDf<-tibble(date=as_date(timeBnd[1,],origin="1850-01-01"))%>%
      mutate(year=year(date))

    for (i in 1:nrow(datesDf)){ # iterate over years

      tempDf <- expand.grid(lat=lat[1,],lon=lon[,1])%>%
        mutate(lon=ifelse(lon>180,-(360-lon),lon),
               o2=as.vector(t(datVar[,,3,i])))%>% # 3 is depth of 25m
        mutate(year=as.numeric(datesDf[i,"year"]))

      o2Dfshallow = o2Dfshallow %>% bind_rows(tempDf)
    }

    ## Deeper: >425m

    o2DfDeep<-NULL

    for (i in 1:nrow(datesDf)){ # iterate over years

      tempDf <- expand.grid(lat=lat[1,],lon=lon[,1])%>%
        mutate(lon=ifelse(lon>180,-(360-lon),lon),
               o2=as.vector(t(datVar[,,26,i])))%>% # 26 is depth of 427.3156m
        mutate(year=as.numeric(datesDf[i,"year"]))

      o2DfDeep = o2DfDeep %>% bind_rows(tempDf)
    }
  }

  ### calculate changes to 2060 and 2100
  if(var %in% c("tos","sos")) {

    varDiffs<-aveDf%>%
      filter(year==2020)%>%
      rename(var2020=aveVar)%>%
      bind_cols(aveDf%>%
                  filter(year==2060)%>%
                  select(var2060=aveVar),
                aveDf%>%
                  filter(year==2100)%>%
                  select(var2100=aveVar))%>%
      mutate(diff2060=var2060-var2020,
             diff2100=var2100-var2020)
  } else {

    o2DiffsSh<-o2Dfshallow%>%
      filter(year==2020)%>%
      rename(o22020=o2)%>%
      bind_cols(o2Dfshallow%>%
                  filter(year==2060)%>%
                  select(o22060=o2),
                o2Dfshallow%>%
                  filter(year==2100)%>%
                  select(o22100=o2))%>%
      mutate(diff2060=o22060-o22020,
             diff2100=o22100-o22020)

    o2DiffsDp<-o2DfDeep%>%
      filter(year==2020)%>%
      rename(o22020=o2)%>%
      bind_cols(o2DfDeep%>%
                  filter(year==2060)%>%
                  select(o22060=o2),
                o2DfDeep%>%
                  filter(year==2100)%>%
                  select(o22100=o2))%>%
      mutate(diff2060=o22060-o22020,
             diff2100=o22100-o22020)

  }

  ### set up plot area and shp files

  if(is.na(corr)) {

    mpa<-mpaShp%>%
      st_set_crs(4326)

    siteBB<-st_bbox(mpa)

    ploArea<-st_bbox(c(xmin=as.numeric(siteBB$xmin-buff),
                       xmax=as.numeric(siteBB$xmax+buff),
                       ymin=as.numeric(siteBB$ymin-buff),
                       ymax=as.numeric(siteBB$ymax+buff)),crs = st_crs(mpa))
  } else {

    mpa<-mpaShp%>%
      st_cast(to="POINT")%>% # convert polygons to points
      mutate(lon=st_coordinates(geometry)[,"X"] ,
             lat=st_coordinates(geometry)[,"Y"])%>%
      mutate(lon2=case_when(lon < 0 ~ 180+(180+lon),
                            TRUE ~ lon))%>%
      st_drop_geometry()%>%
      sfheaders::sf_multipolygon(x="lon2",y="lat",multipolygon_id=colnames(mpaShp)[1])%>%
      st_set_crs(4326)

    siteBB<-st_bbox(mpa)

    ploArea<-st_bbox(c(xmin=as.numeric(ifelse(siteBB$xmin<0,180+(180+siteBB$xmin),siteBB$xmin)-buff),
                       xmax=as.numeric(ifelse(siteBB$xmax<0,180+(180+siteBB$xmax),siteBB$xmax)+buff),
                       ymin=as.numeric(siteBB$ymin-buff),
                       ymax=as.numeric(siteBB$ymax+buff)),crs = st_crs(mpa))
  }

  ### Clip var data to site

  if(is.na(corr)) {
    if(var %in% c("tos","sos")){

      varSite<-varDiffs%>%
        filter(lon>as.numeric(ploArea$xmin),
               lon<as.numeric(ploArea$xmax),
               lat>as.numeric(ploArea$ymin),
               lat<as.numeric(ploArea$ymax))
    } else {
      o2SiteSh<-o2DiffsSh%>%
        filter(lon>as.numeric(ploArea$xmin),
               lon<as.numeric(ploArea$xmax),
               lat>as.numeric(ploArea$ymin),
               lat<as.numeric(ploArea$ymax))%>%
        mutate(ave2020=mean(o22020,na.rm=T),
               o2perc60=(diff2060/ave2020)*100)

      o2SiteDp<-o2DiffsDp%>%
        filter(lon>as.numeric(ploArea$xmin),
               lon<as.numeric(ploArea$xmax),
               lat>as.numeric(ploArea$ymin),
               lat<as.numeric(ploArea$ymax))%>%
        mutate(ave2020=mean(o22020,na.rm=T),
               o2perc60=(diff2060/ave2020)*100)
    }
  } else {

    if(var %in% c("tos","sos")){

      varSite<-varDiffs%>%
        mutate(lon=case_when(lon < 0 ~ 180+(180+lon),
                             TRUE ~ lon))%>%
        filter(lon>as.numeric(ploArea$xmin),
               lon<as.numeric(ploArea$xmax),
               lat>as.numeric(ploArea$ymin),
               lat<as.numeric(ploArea$ymax))
    } else {
      o2SiteSh<-o2DiffsSh%>%
        mutate(lon=case_when(lon < 0 ~ 180+(180+lon),
                             TRUE ~ lon))%>%
        filter(lon>as.numeric(ploArea$xmin),
               lon<as.numeric(ploArea$xmax),
               lat>as.numeric(ploArea$ymin),
               lat<as.numeric(ploArea$ymax))%>%
        mutate(ave2020=mean(o22020,na.rm=T),
               o2perc60=(diff2060/ave2020)*100)

      o2SiteDp<-o2DiffsDp%>%
        mutate(lon=case_when(lon < 0 ~ 180+(180+lon),
                             TRUE ~ lon))%>%
        filter(lon>as.numeric(ploArea$xmin),
               lon<as.numeric(ploArea$xmax),
               lat>as.numeric(ploArea$ymin),
               lat<as.numeric(ploArea$ymax))%>%
        mutate(ave2020=mean(o22020,na.rm=T),
               o2perc60=(diff2060/ave2020)*100)
    }
  }

  ### crop base layers from world map

  if(is.na(corr)) {
    geoShp<-map_data("world")%>%
      st_as_sf(.,
               coords=c("long","lat"),
               crs=4326)%>%#"+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")%>%
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
      st_make_valid()%>%
      st_set_crs(4326)
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

  ### plot data
  if(var %in% c("tos","sos")) {

    if(any(varSite$diff2060 > 0) & any(varSite < 0)) {
      sclFllGrd<-scale_fill_gradient2(low = "#0a6165",high = "#da4325",mid = "#e7e1bc",midpoint = 0,
                                      name=ifelse(var=="tos","degC","ppt"))
    } else {

      if(all(varSite$diff2060>=0)) {
        sclFllGrd<-scale_fill_gradient(high = "#da4325",low = "#e7e1bc",
                                       name=ifelse(var=="tos","degC","ppt"))
      } else {
        sclFllGrd<-scale_fill_gradient(low = "#0a6165",high = "#e7e1bc",
                                       name=ifelse(var=="tos","degC","ppt"))
      }
    }

    sitePlt<-ggplot()+
      geom_tile(data=varSite,aes(x=lon,y=lat,fill=diff2060),width=1,height=1)+
      geom_sf(data=mpa,color="grey21",fill=NA)+
      geom_sf(data=geoShp)+
      # geom_sf(data=mrRev)+
      sclFllGrd+
      scale_x_continuous(expand = c(0,0))+
      scale_y_continuous(expand = c(0,0),limits = c(ploArea$ymin,ploArea$ymax+0.1))+
      theme(legend.direction = "vertical",legend.position = "right")+
      ggtitle(paste("Projected Change in Sea Surface",ifelse(var=="tos","Temperature","Salinity")),subtitle = "2020 to 2060")+
      theme_bw()

    return(sitePlt)
  } else {

    if(any(o2SiteSh$diff2060 > 0) & any(o2SiteSh < 0)) {
      sclFllGrdSh<-scale_fill_gradient2(low = "#0a6165",high = "#da4325",mid = "#e7e1bc",midpoint = 0,
                                        name="% change in\n[oxygen]")
    } else {

      if(all(o2SiteSh$diff2060>=0)) {
        sclFllGrdSh<-scale_fill_gradient(high = "#da4325",low = "#e7e1bc",
                                         name="% change in\n[oxygen]")
      } else {
        sclFllGrdSh<-scale_fill_gradient(low = "#0a6165",high = "#e7e1bc",
                                         name="% change in\n[oxygen]")
      }
    }

    o2DiffSh<-ggplot()+
      geom_tile(data=o2SiteSh,aes(x=lon,y=lat,fill=o2perc60),width=1,height=1)+
      geom_sf(data=mpa,color="white",fill=NA)+
      geom_sf(data=geoShp)+
      # geom_sf(data=mrRev)+
      sclFllGrdSh+
      scale_x_continuous(expand = c(0,0))+
      scale_y_continuous(expand = c(0,0),limits = c(ploArea$ymin,ploArea$ymax+0.1))+
      theme(legend.direction = "vertical",legend.position = "right")+
      ggtitle("Projected Change in Oxygen Concentration: shallow (25m)",subtitle = "2020 to 2060")+
      theme_bw()

    if(any(o2SiteDp$diff2060 > 0) & any(o2SiteDp < 0)) {
      sclFllGrdDp<-scale_fill_gradient2(low = "#0a6165",high = "#da4325",mid = "#e7e1bc",midpoint = 0,
                                        name="% change in\n[oxygen]")
    } else {

      if(all(o2SiteDp$diff2060>=0)) {
        sclFllGrdDp<-scale_fill_gradient(high = "#da4325",low = "#e7e1bc",
                                         name="% change in\n[oxygen]")
      } else {
        sclFllGrdDp<-scale_fill_gradient(low = "#0a6165",high = "#e7e1bc",
                                         name="% change in\n[oxygen]")
      }
    }

    o2DiffDp<-ggplot()+
      geom_tile(data=o2SiteDp,aes(x=lon,y=lat,fill=o2perc60),width=1,height=1)+
      geom_sf(data=mpa,color="white",fill=NA)+
      geom_sf(data=geoShp)+
      # geom_sf(data=mrRev)+
      sclFllGrdDp+
      scale_x_continuous(expand = c(0,0))+
      scale_y_continuous(expand = c(0,0),limits = c(ploArea$ymin,ploArea$ymax+0.1))+
      theme(legend.direction = "vertical",legend.position = "right")+
      ggtitle("Projected Change in Oxygen Concentration: Deep (>425m)",subtitle = "2020 to 2060")+
      theme_bw()

    return(list(o2DiffSh,o2DiffDp))
  }

}
