#' Title
#'
#' @param mpaShp shp file of the existing or proposed MPA
#' @param buff distance around mpa shape file for which analysis is conducted (to calculate extent box, "ploArea"). In decimal degrees (numeric)
#' @param siteName character string of site name (no spaces) for file name creation
#' @param nowFut one of "now" or "fut" to indicate which dataset to pull from: now for current ranges, fut for 2050 projections
#' @param corr Indicates whether longitude "lonFx" needs to be corrected for crossing the 180 meridian, if no use NA
#' @param datFolder file path to where aqm data are stored
#'
#' @return writes a spatial features dataframe with aquamaps data from the indicated dataset (now/fut) for the defined plot area
#' @export
#'
#' @examples
#'
#'aqmClip(mpaShp=revMPA,buff=10,siteName="rev",datFolder="../../../ccAssmts/data")
aqmClip<-function(mpaShp,
                  buff,
                  siteName,
                  nowFut=c("now","fut"),
                  corr=c(NA,"lonFx"),
                  datFolder) {

  mpa<-mpaShp%>%
    st_set_crs("EPSG 4326")#"+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

  ### indicate which data file
  if(nowFut=="now"){
    aqData<-paste(datFolder,"hcaf_species_native-001.csv",sep="/")
  } else {
    aqData<-paste(datFolder,"rcp85_hcaf_species_native_2050.csv",sep="/")
  }

  ### download the data
  # datFile<-"data/tempDat/aquamapsDat.csv" #where to put it
  #
  # drive_download(aqData, #csv file
  #                path=datFile,
  #                overwrite = TRUE)
  #
  # aqDat<-read_csv(datFile,
  #                 col_select = c(SpeciesID,CenterLat,CenterLong,Probability))
  #
  # ### delete large file once the object has been created
  #
  # if (file.exists(datFile)) {
  #   unlink(datFile)
  #   cat(paste("The file has been deleted: ",datFile))
  # }

  ## assuming I'm working on the server: direct pull in of the data

  aqDat<-read_csv(aqData)

  ### establish site coordinates
  siteBB<-st_bbox(mpaShp)

  if(is.na(corr)) {
    ploArea<-st_bbox(c(xmin=as.numeric(siteBB$xmin-buff),
                       xmax=as.numeric(siteBB$xmax+buff),
                       ymin=as.numeric(siteBB$ymin-buff),
                       ymax=as.numeric(siteBB$ymax+buff)),crs = st_crs(mpa))
  } else {
    ploArea<-st_bbox(c(xmin=as.numeric(ifelse(siteBB$xmin<0,180+(180+siteBB$xmin),siteBB$xmin)-buff),
                       xmax=as.numeric(ifelse(siteBB$xmax<0,180+(180+siteBB$xmax),siteBB$xmax)+buff),
                       ymin=as.numeric(siteBB$ymin-buff),
                       ymax=as.numeric(siteBB$ymax+buff)),crs = st_crs(mpa))
  }

  ### Clip aquamaps data to site

  if(is.na(corr)) {

    aqSiteDat<-aqDat%>%
      filter(CenterLong>as.numeric(ploArea$xmin),
             CenterLong<as.numeric(ploArea$xmax),
             CenterLat>as.numeric(ploArea$ymin),
             CenterLat<as.numeric(ploArea$ymax))%>%
      st_as_sf(coords=c("CenterLong","CenterLat"),
               crs=4326)

  } else {

    aqSiteDat<-aqDat%>%
      mutate(lon=case_when(CenterLong < 0 ~ 180+(180+CenterLong),
                           TRUE ~ CenterLong))%>%
      rename(lat=CenterLat)%>%
      filter(lon>ploArea$xmin & lon<ploArea$xmax,
             lat>ploArea$ymin & lat<ploArea$ymax)%>%
      st_as_sf(coords=c("lon","lat"),
               crs=4326)

  }

  fileNm<-paste("data/",siteName,"Aq",nowFut,".shp",sep="")

  ### write to new file
  st_write(aqSiteDat,fileNm)

  # print(fileNm)
  # return(aqSiteDat)

}
