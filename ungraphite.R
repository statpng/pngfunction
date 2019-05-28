ungraphite = function(species="hsapiens", to="SYMBOL", out=c("pathway","eframe","adjm"), edge = NULL){

  library(graphite)
  
  require(graphite)
  require(igraph) ## adjacency matrix
  
  
  out = match.arg(out)
  species = match.arg(species)
  to = match.arg(to)
  
  db <- graphite::pathwayDatabases()
  db <- as.character( db[ db[,1] == species , 2] )
  
  Version = packageVersion("graphite")
  
  
  if( paste0("hsapiens_",to,"_",Version,".RData") %in% list.files() ){
    load(paste0("hsapiens_",to,".RData"))
    
  } else {
    hsapiens <- lapply( db , function(x) graphite::pathways( species , x ) )
    cat("Convert Identifier \n")
    hsapiens_to <- lapply( hsapiens , function(x) convertIdentifiers(x=x , to = "SYMBOL") )
    names(hsapiens_to) = db
    save( hsapiens_to , file=paste0("hsapiens_",to,"_",Version,".RData") ) 
  }
  
  cat("Creating edge frame \n")
  if(out=="pathway") return(hsapiens_to)
  
  if( paste0("edge_hsapiens_",to,"_",Version,".RData") %in% list.files() ){
    
    load( paste0("edge_hsapiens_",to,"_",Version,".RData") )
    
  } else {
    
    eframe <- lapply( hsapiens_to, function(x) lapply(x, function(y) y@edges[,1:2]) )
    
    eframe <- lapply( eframe , function(x) do.call( rbind , x ) )
    
    eframe <- do.call( rbind , eframe )
    rownames(eframe) = NULL
    
    eframe <- eframe[!duplicated(eframe),]
    
    eframe <- eframe[ !apply( eframe , 1 , function(x) x[1] == x[2] ) , ]
    
    eframe <- t( apply( eframe , 1 , sort ) )
    
    eframe <- eframe[!duplicated(eframe),]
    edge_hsapiens <- as.data.frame(eframe)
    colnames(edge_hsapiens) <- c("source","dest")
    
    save( edge_hsapiens , file=paste0("edge_hsapiens_",to,"_",Version,".RData") )
    
  }
  
  if(out=="eframe") return(edge_hsapiens)
  
  if( out=="adjm" & is.null(edge) ){
    adj <- as_adj( graph_from_data_frame(edge_hsapiens, directed = FALSE, vertices = NULL), type="both" )
  } else if(out=="adjm" & !is.null(edge) ){
    
    eframe <- edge_hsapiens[ apply( apply( edge_hsapiens , 2 , function(x) x %in% edge ) , 1 , all ) , ]
    
    adj <- as_adj( graph_from_data_frame(eframe, directed = FALSE, vertices = NULL), type="both" )
    
    #adj <- adjacencyList2Matrix( eframe , square=TRUE)
  }
  
  res = list( eframe=eframe, adjm = adj )
  
  if(out=="adjm") return(res)
  
}

# pathway = ungraphite(out="pathway") ## Human(homo sapiens) Gene이 저장된 DB 6곳(KEGG, Reactome 등)의 graphite 클래스를 불러옴
# eframe = ungraphite(out="eframe" ) ## 6곳의 DB 자료에서 Edge만 추출하여 방향성(<x1,x2> == <x2,x1>), 중복(Duplicated),  등을 제거
# adjm = ungraphite(out="adjm", edge = c("ABCB4","ABTB2","ACO1","AKAP13","ANXA11","C1R","PDE11A","PDE4C","PDE4D","COL18A1","MMP3") )
## 특정 nodes에 대한 데이터만 추출.(argument "edge"가 사실은 "nodes"로 수정해야함)

# 각 DB에 저장된 Gene 형식이 다름; "symbol","entrez ID" 등(graphite::convertIdentifiers 함수, Annotation packge 참고)
# 위 형식을 통일해주기 위한 함수(convertIdentifiers)에서 시간이 오래 걸림(약 5분).
# 따라서 함수를 처음 실행하면 작업디렉토리에 작업결과를 저장한 RData를 생성.
