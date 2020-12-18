# Library Load
rm(list=ls())

# 종합기상관측 (Automated Synoptic Observing System; ASOS)
# https://data.kma.go.kr/data/grnd/selectAsosRltmList.do?pgmNo=36
# How to get a data from open API
# https://www.2e.co.kr/news/articleView.html?idxno=210202

api_key <- "vSOirPtpYtebEsKbUTv5Wyy2OZ1Kl453aYa6FzKy5UUN0BrKQzP%2BS7HFbvJkbSomRN5yZMacxzj6VhteDtMwDg%3D%3D"
serviceURL <- "http://apis.data.go.kr/1360000/AsosDalyInfoService/getWthrDataList"
numOfRows <- 100
dataType <- "XML"
dataCd <- "ASOS"
dateCd <- "DAY"
startDt <- 19810101
endDt <- 20191231
stnIds <- 108

png.get_API <- function(api_key,
                        serviceURL,
                        numOfRows,
                        dataType,
                        dataCd,
                        dateCd,
                        startDt,
                        endDt,
                        stnIds){
  
  url <- paste0( serviceURL, 
                 "?serviceKey=", api_key,
                 "&pageNo=", 1,
                 "&numOfRows=", numOfRows,
                 "&dataType=", dataType,
                 "&dataCd=", dataCd,
                 "&dateCd=", dateCd,
                 "&startDt=", startDt,
                 "&endDt=", endDt,
                 "&stnIds=", stnIds)
  
  xmlDoc <- xmlTreeParse(url, useInternalNodes = TRUE, encoding = "UTF-8")
  rootNode <- xmlRoot(xmlDoc)
  numOfRows <- as.numeric(xpathSApply(rootNode, "//numOfRows", xmlValue))
  totalCount <- as.numeric(xpathSApply(rootNode, "//totalCount", xmlValue))
  loopCount <- round( totalCount/numOfRows, 0 )
  if( loopCount * numOfRows < totalCount ){
    loopCount <- loopCount + 1
  }
  
  totalData <- NULL
  for(i in 1:loopCount){
    url <- paste0( serviceURL, 
                   "?serviceKey=", api_key,
                   "&pageNo=", 1,
                   "&numOfRows=", numOfRows,
                   "&dataType=", dataType,
                   "&dataCd=", dataCd,
                   "&dateCd=", dateCd,
                   "&startDt=", startDt,
                   "&endDt=", endDt,
                   "&stnIds=", stnIds)
    
    xmlDoc <- xmlTreeParse(url, useInternalNodes = TRUE, encoding = "UTF-8")
    rootNode <- xmlRoot(xmlDoc)
    
    xmlData <- xmlToDataFrame(nodes = getNodeSet(rootNode, "//item"))
    
    totalData <- rbind( totalData, xmlData)
  }
  
}
