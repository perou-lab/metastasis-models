# Functions

readarray<-function(dataFile,designFile=NA,hr=1,impute=T,method="mean"){
  
  headerRows <- hr
  
  x<-read.table(dataFile,sep="\t",header=F,fill=T,stringsAsFactors=FALSE)
  
  if(headerRows==1){
    sampleNames<-as.vector(t(x[1,-1]))
    x<-x[-1,]
    classes<-NULL
    ids<-x[,1]
    xd<-x[,-1]
    xd<-apply(xd,2,as.numeric)
    xd<-collapseIDs(xd,ids,method)	
  }else{
    sampleNames<-as.vector(t(x[1,-1]))
    x<-x[-1,]
    
    classes<-x[1:(headerRows-1),]
    dimnames(classes)[[1]]<-classes[,1]
    classes<-classes[,-1]
    classes[classes==""]<-NA
    classes<-t(classes)
    rownames(classes)<-sampleNames
    classes<-as.data.frame(classes)
    
    xd<-x[(-1:-(headerRows-1)),]
    ids<-as.vector(t(xd[,1]))
    xd<-xd[,-1]
    xd<-apply(xd,2,as.numeric)
    xd<-collapseIDs(xd,ids,method)
  }
  
  features<- dim(xd)[1]
  samples<- dim(xd)[2]
  geneNames<-rownames(xd)
  xd<-apply(xd,2,as.numeric)
  rownames(xd)<-geneNames
  colnames(xd)<-sampleNames
  
  if(!is.na(designFile)){
    x<-read.table(designFile,sep="\t",header=T,row.names=1,fill=T,,stringsAsFactors=FALSE)
    xd<-xd[,sort.list(colnames(xd))]
    xd<-xd[,colnames(xd) %in% rownames(x)]
    x<-x[rownames(x) %in% colnames(xd),]
    x<-x[sort.list(rownames(x)),]
    classes<-as.data.frame(x)
  }
  
  if(sum(apply(xd,2,is.na))>0 & impute){
    library(impute)
    allAnn<-dimnames(xd)
    data.imputed<-impute.knn(as.matrix(xd))$data
    xd<-data.imputed[1:features,]
    dimnames(xd)<-allAnn
  }
  
  return(list(xd=xd, classes=classes, nfeatures=features, nsamples=samples, fnames=geneNames, snames=sampleNames))
}


collapseIDs<-function(x,allids=row.names(x),method="mean"){
  
  allids<-as.vector(allids)
  x.col<- NULL
  
  ids<-as.character(unique(allids[duplicated(allids)]))
  
  if(length(ids)== 0){ 
    dimnames(x)[[1]]<-allids
    return(x) 
  }
  
  
  for(i in 1:length(ids)){
    indices <- allids==ids[i] 
    if(method=="mean"){
      vals<-apply(x[indices,],2,mean,na.rm=T)
    }
    if(method=="median"){
      vals<-apply(x[indices,],2,median,na.rm=T)
    }
    if(method=="max"){
      vals<-apply(x[indices,],2,max,na.rm=T)
    }
    if(method=="sum"){
      vals<-apply(x[indices,],2,sum,na.rm=T)
    }
    if(method=="stdev"){   
      temp<- x[indices,]
      stdevs<- apply(temp,1,sd,na.rm=T)
      vals<- temp[match(max(stdevs),stdevs),]
    }
    if(method=="iqr"){   
      temp<- x[indices,]
      iqrs<- apply(temp,1,function(x){quantile(x,.75,na.rm=T)-quantile(x,.25,na.rm=T)})
      vals<- temp[match(max(iqrs),iqrs),]
    }
    x.col <- rbind(x.col,vals)
  }
  
  
  non.overlapping.ids <- !(allids %in% ids) 
  
  x.col <- rbind(x.col,x[non.overlapping.ids,])
  dimnames(x.col)<- list(c(ids,allids[non.overlapping.ids]),dimnames(x)[[2]])
  
  return(x.col)
  
}