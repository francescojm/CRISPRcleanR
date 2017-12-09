## exported, documented, vignetted functions
#### Analysis
ccr.NormfoldChanges<-function(filename,
                              display=TRUE,
                              saveToFig=FALSE,
                              outdir='./',
                              min_reads=30,
                              EXPname='',
                              libraryAnnotation,
                              ncontrols=1){
    
    counts<-read.table(filename,sep='\t',header=TRUE,stringsAsFactors = FALSE)
    counts<-counts[is.element(counts$sgRNA,rownames(libraryAnnotation)),]
    
    if(saveToFig){
        display<-TRUE
        pdf(paste(outdir,EXPname,'_normCounts.pdf',sep=''),width = 10,height = 10)
    }
    
    if(display){
        par(mfrow=c(1,2))
        ccr.boxplot(counts[,3:ncol(counts)],
                    main = paste(EXPname,'Raw sgRNA counts'),
                    names = c(paste('CTRL',1:ncontrols),paste('library r',1:(ncol(counts)-2-ncontrols))))
    }
    
    numd<-
        counts[is.element(counts$sgRNA,rownames(libraryAnnotation)),3:ncol(counts)]
    
    if (ncontrols==1){
        IDX<-which(numd[,1]>=min_reads)    
    }else{
        IDX<-which(rowMeans(numd[,1:ncontrols])>=min_reads)
    }
    
    numd<-numd[IDX,]
    counts<-counts[IDX,]
    
    normFact<-t(matrix(rep(colSums(numd),nrow(numd)),ncol(counts)-2,nrow(numd)))
    
    numd<-numd/normFact*10000000
    
    normed<-cbind(counts[,1:2],numd)
    
    if(display){
        ccr.boxplot(numd,main = paste(EXPname,'normalised sgRNA counts'),
                    names = c(paste('CTRL',1:ncontrols),paste('library r',1:(ncol(counts)-2-ncontrols))))
    }
    
    if(saveToFig){
        dev.off()
    }
    
    
    nsamples<-ncol(counts)-(2+ncontrols)
    
    for ( i in 1:nsamples){
        
        if(ncontrols==1){
            c_foldchanges<-log2((normed[,3+i]+0.5)/(normed[,3]+0.5))    
        }else{
            c_foldchanges<-log2((normed[,2+ncontrols+i]+0.5)/(rowMeans(normed[,3:(2+ncontrols)])+0.5))
        }
        
        c_foldchanges<-matrix(c_foldchanges,length(c_foldchanges),1,
                              dimnames = list(normed$sgRNA,colnames(normed)[2+ncontrols+i]))
        
        if (i == 1){
            foldchanges<-c_foldchanges
        }else{
            foldchanges<-cbind(foldchanges,c_foldchanges )
        }
    }
    
    if(saveToFig){
        pdf(paste(outdir,EXPname,'_fcs.pdf',sep=''),width = 8,height = 10)    
    }
    
    if(display){
        par(mfrow=c(1,1))
        ccr.boxplot(foldchanges,
                    main = paste(EXPname,'log Fold Changes'),
                    names = paste('library r',1:(ncol(counts)-2-ncontrols)))
        if(saveToFig){
            dev.off()    
        }
        
    }
    
    foldchanges<-cbind(normed[,1:2],foldchanges)
    
    save(normed,file=paste(outdir,EXPname,'_normCounts.Rdata',sep=''))
    save(foldchanges,file=paste(outdir,EXPname,'_foldChanges.Rdata',sep=''))
    
    return(list(norm_counts=normed,logFCs=foldchanges))
    
}

ccr.logFCs2chromPos<-function(foldchanges,libraryAnnotation){
    
    sgRNAsIds<-foldchanges$sgRNA
    
    genes<-
        as.character(libraryAnnotation[sgRNAsIds,"GENES"])
    
    chrN<-
        as.character(libraryAnnotation[sgRNAsIds,"CHRM"])
    
    chrN[which(chrN=='X')]<-'23'
    chrN[which(chrN=='Y')]<-'24'
    chrN<-as.numeric(chrN)
    
    startp<-as.numeric(as.character(libraryAnnotation[sgRNAsIds,"STARTpos"]))
    endp<-as.numeric(as.character(libraryAnnotation[sgRNAsIds,"ENDpos"]))
    
    position<-round(startp+(endp-startp)/2)
    
    if (ncol(foldchanges)>3){
        converted<-
            data.frame(chrN,startp,endp,genes,rowMeans(foldchanges[,3:ncol(foldchanges)]),stringsAsFactors = FALSE)
    }else{
        converted<-
            data.frame(chrN,startp,endp,genes,foldchanges[,3:ncol(foldchanges)],stringsAsFactors = FALSE)
    }
    
    rownames(converted)<-foldchanges$sgRNA
    
    converted<-converted[order(converted$chrN,converted$startp),]
    
    colnames(converted)[5]<-'avgFC'
    
    colnames(converted)[1]<-'CHR'
    BP<-converted$startp+(converted$endp-converted$startp)/2
    
    converted<-cbind(converted,BP)
    
    return(converted)
}
ccr.cleanChrm<-function(gwSortedFCs,
                        CHR,
                        display=TRUE,
                        label='',
                        saveTO=NULL,
                        min.ngenes=3,
                        ignoredGenes=NULL){
    
    ID<-which(gwSortedFCs$CHR==CHR)
    gwSortedFCs<-gwSortedFCs[ID,]
    
    my.CNA.object <- CNA(cbind(gwSortedFCs$avgFC),
                         gwSortedFCs$CHR,
                         gwSortedFCs$BP,
                         data.type="logratio",sampleid=paste(label,'Chr',CHR,'sgRNA FCs'))
    
    
    
    my.smoothed.CNA.object <- smooth.CNA(my.CNA.object)
    my.segment.smoothed.CNA.object <- segment(my.smoothed.CNA.object, verbose=1)
    
    nsegments<-
        nrow(my.segment.smoothed.CNA.object$output)
    
    nGeneInSeg<-vector()
    guides<-vector()
    newFC<-gwSortedFCs$avgFC
    correction<-rep(0,length(newFC))
    
    for (i in 1:nsegments){
        idxs<-my.segment.smoothed.CNA.object$segRows[i,1]:my.segment.smoothed.CNA.object$segRows[i,2]
        includedGenes<-unique(gwSortedFCs[idxs,'genes'])
        includedGuides<-ID[range(idxs)]
        
        if(length(ignoredGenes)>0){
            nGeneInSeg[i]<-length(setdiff(includedGenes,ignoredGenes))    
        }else{
            nGeneInSeg[i]<-length(includedGenes)
        }
        
        
        if(nGeneInSeg[i]>=min.ngenes){
            
            orSign<-sign(newFC[idxs])
            newFC[idxs]<-newFC[idxs]-mean(newFC[idxs])
            
            correction[idxs]<- -sign(mean(newFC[idxs]))
        }
        guides[i]<-paste(includedGuides,collapse=', ')
    }
    
    norm.my.CNA.object <- CNA(cbind(newFC),
                              gwSortedFCs$CHR,
                              gwSortedFCs$BP,
                              data.type="logratio",sampleid=paste(label,'Chr',CHR,'sgRNA FCs - post CRISPRcleanR'))
    
    norm.my.smoothed.CNA.object <- smooth.CNA(norm.my.CNA.object)
    norm.my.segment.smoothed.CNA.object <- segment(norm.my.smoothed.CNA.object, verbose=1)
    
    if(!display){
        saveTO<-NULL
    }
    
    if(length(saveTO)){
        display<-TRUE
        path<-paste(saveTO,label,'/',sep='')
        if(!file.exists(path)){
            dir.create(path)
        }
        pdf(paste(path,CHR,'.pdf',sep=''),width = 7.5,height = 7.5)
    }
    
    if(display){
        par(mfrow=c(2,1))
        plot(my.segment.smoothed.CNA.object,
             main=paste(label,'Chr',CHR,'sgRNA FCs'),
             segcol='black',xmaploc=FALSE,xlab='',ylab='logFC')
        
        plot(norm.my.segment.smoothed.CNA.object,main= paste(label,'Chr',CHR,'sgRNA FCs, post CRISPRcleanR'),
             segcol='black',xmaploc=FALSE,xlab='sgRNA index',ylab='logFC')    
    }
    
    if(length(saveTO)){
        dev.off()    
    }
    
    
    correctedFC<-newFC
    gwSortedFCs<-cbind(gwSortedFCs,correction,correctedFC)
    
    regions<-my.segment.smoothed.CNA.object
    regions<-regions$output
    
    regions<-cbind(regions[,2:ncol(regions)],guides)
    
    colnames(regions)<-c('CHR','startp','endp','n.sgRNAs','avg.logFC','guideIdx')
    
    res<-list(correctedFCs=gwSortedFCs,regions=regions)
    return(res)
    
}

ccr.GWclean<-function(gwSortedFCs,label='',display=TRUE,saveTO=NULL,ignoredGenes=NULL,min.ngenes=3){
    
    CHRs<-as.character(sort(unique(gwSortedFCs$CHR)))
    
    for (i in 1:length(CHRs)){
        CHR<-CHRs[i]
        res<-ccr.cleanChrm(gwSortedFCs = gwSortedFCs,
                           CHR = CHR,
                           display = display,
                           label = label,
                           saveTO=saveTO,
                           ignoredGenes = ignoredGenes,
                           min.ngenes=min.ngenes)
        
        if(i == 1){
            corrected_logFCs<-res$correctedFCs
            segments<-res$regions
        }else{
            corrected_logFCs<-rbind(corrected_logFCs,res$correctedFCs)
            segments<-rbind(segments,res$regions)
        }
    }
    
    return(list(corrected_logFCs=corrected_logFCs,segments=segments,SORTED_sgRNAs=rownames(gwSortedFCs)))
    
}
ccr.correctCounts<-function(CL,normalised_counts,correctedFCs_and_segments,
                            libraryAnnotation,
                            minTargetedGenes=3,
                            OutDir='./',
                            ncontrols=1){
    
    normalised_counts<-normalised_counts[which(!is.na(normalised_counts[,1])),]
    rownames(normalised_counts)<-normalised_counts$sgRNA
    
    numdata<-normalised_counts[,3:ncol(normalised_counts)]
    rownames(numdata)<-normalised_counts$sgRNA
    
    segments<-correctedFCs_and_segments$segments
    
    segment_guides<-list()
    
    nsegments<-nrow(segments)
    
    correctedCounts<-numdata
    
    for (i in 1:nsegments){
        print(paste(CL,': correcting counts on segment ',i,' (of ',nsegments,')',sep=''))
        current_idxs<-as.character(segments$guideIdx[i])
        current_idxs<-as.numeric(unlist(str_split(current_idxs,', ')[[1]]))
        segment_guides[[i]]<-correctedFCs_and_segments$SORTED_sgRNAs[current_idxs[[1]]:current_idxs[[2]]]
        
        ntarg<-length(unique(as.character(libraryAnnotation[segment_guides[[i]],'GENES'])))
        
        guides<-segment_guides[[i]]
        if (ntarg>=minTargetedGenes){
            FC<-correctedFCs_and_segments$corrected_logFCs[guides,'avgFC']
            
            if(ncontrols==1){
                c<-numdata[unlist(guides),1]    
            }else{
                c<-rowMeans(numdata[unlist(guides),1:ncontrols])
            }
            
            
            N<-correctedFCs_and_segments$corrected_logFCs[guides,'correctedFC']
            Cf<-N-FC
            
            #reverted<-(N^2)+(2*N*log2(c))+c
            reverted<-c*2^N
            
            nreps<-ncol(numdata)-ncontrols
            
            if(nreps>1){
                proportions<-(numdata[guides,(ncontrols+1):ncol(numdata)]+1)/
                    (matrix(rep(rowSums(numdata[guides,(ncontrols+1):ncol(numdata)]+1),
                                nreps),length(guides),nreps))
                revertedCounts<-(reverted*nreps)*proportions
            }else{
                revertedCounts<-reverted
            }
            
            
            correctedCounts[guides,(ncontrols+1):ncol(numdata)]<-revertedCounts
        }
    }
    
    adjusted<-as.data.frame(cbind(normalised_counts$sgRNA,normalised_counts$gene,correctedCounts),stringAsFactors=FALSE)
    colnames(adjusted)<-colnames(normalised_counts)
    rownames(adjusted)<-NULL
    
    correctedCounts<-adjusted
    
    save(correctedCounts,file=paste(OutDir,CL,'_correctedCounts.RData',sep=''))
    return(correctedCounts)
}

#### Utils
ccr.genes2sgRNAs<-function(libraryAnnotation,genes){
    notIncludedGenes<-genes[which(!is.element(genes,libraryAnnotation$GENES))]
    if(length(notIncludedGenes)>0){
        warning(paste('No sgRNAs targeting the following genes in this library:',paste(notIncludedGenes,collapse=', ')))    
    }
    sgRNAs<-rownames(libraryAnnotation)[which(is.element(libraryAnnotation$GENES,genes))]
    return(sgRNAs)
}
ccr.geneMeanFCs<-function(sgRNA_FCprofile,libraryAnnotation){
    FCsprofile<-sgRNA_FCprofile
    FCsprofile<-aggregate(sgRNA_FCprofile~libraryAnnotation[names(sgRNA_FCprofile),'GENES'], FUN=mean)
    nn<-as.character(FCsprofile$`libraryAnnotation[names(sgRNA_FCprofile), \"GENES\"]`)
    FCsprofile<-FCsprofile$sgRNA_FCprofile
    names(FCsprofile)<-as.character(nn)
    return(FCsprofile)
}
ccr.get.gdsc1000.AMPgenes<-function(cellLine,
                               minCN=8,
                               exact=FALSE,
                               GDSC.geneLevCNA=NULL){
    
    if(is.null(GDSC.geneLevCNA)){
        data(GDSC.geneLevCNA,envir = environment())    
    }
    
    data(GDSC.CL_annotation)
    
    if(!is.element(cellLine,GDSC.CL_annotation$CL.name) & !is.element(cellLine,GDSC.CL_annotation$COSMIC.ID)){
        stop('Cell line not found: Please provide a valid cell line cosmic identifier or name (see available cell lines here: http://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources//Data/suppData/TableS1E.xlsx)', call. = TRUE, domain = NULL)
    }
    
    if(is.element(cellLine,GDSC.CL_annotation$COSMIC.ID)){
        cid<-cellLine
    }else{
        cid<-as.character(GDSC.CL_annotation$COSMIC.ID[GDSC.CL_annotation$CL.name==cellLine])
    }
    
    
    genes<-rownames(GDSC.geneLevCNA)
    
    if(!exact){
        cnRange<-c(minCN:16)
    }else{
        cnRange<-minCN
    }
    
    if(length(intersect(cid,colnames(GDSC.geneLevCNA)))>0){
        
        amplifiedGenes<-NULL
        
        
        for (ar in cnRange){
            amplifiedGenes<-c(amplifiedGenes,genes[unique(c(grep(paste(',',ar,',',sep=''),
                                                                 GDSC.geneLevCNA[,cid])))])
        }
        
        
        amplifiedGenes<-sort(amplifiedGenes)
        
        GENES<-amplifiedGenes
        CNAval<-GDSC.geneLevCNA[GENES,cid]
        
        CNAval<-unlist(str_split(CNAval,','))
        
        if (length(GENES)>0){
            CN<-CNAval[seq(2,length(GENES)*4,4)]
            
            amplifiedGenes<-cbind(GENES,CN)
            amplifiedGenes<-as.data.frame(amplifiedGenes,stringsAsFactors = FALSE)
            
            colnames(amplifiedGenes)[1]<-'Gene'
            amplifiedGenes[,2]<-as.numeric(amplifiedGenes[,2])
            return(amplifiedGenes)
        }else{
            print('No amplified genes for this cell line according to the selected threshold')
            return(NULL)
        }
    }else{
        print('There is no data available for this cell line in the built in CNA data frame. Provide gene level CN values in input')
        return(NULL)
    }
}
ccr.get.nonExpGenes<-function(cellLine,th=0.05,amplified=FALSE,minCN=8,RNAseq.fpkms=NULL){
    
    data(GDSC.CL_annotation,envir = environment())
    
    if(is.null(RNAseq.fpkms)){
        data(RNAseq.fpkms,envir = environment())   
    }
    
    if(!is.element(cellLine,GDSC.CL_annotation$CL.name) & !is.element(cellLine,GDSC.CL_annotation$COSMIC.ID)){
        stop('Cell line not found: Please provide a valid cell line cosmic identifier or name (see available cell lines here: http://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources//Data/suppData/TableS1E.xlsx)', call. = TRUE, domain = NULL)
    }
    
    if(is.element(cellLine,GDSC.CL_annotation$COSMIC.ID)){
        cid<-cellLine
    }else{
        cid<-as.character(GDSC.CL_annotation$COSMIC.ID[GDSC.CL_annotation$CL.name==cellLine])
    }
    
    
    if(amplified){
        amplifiedGenes<-ccr.get.gdsc1000.AMPgenes(cellLine = cellLine,minCN=minCN)
        amplifiedGenes<-
            unique(amplifiedGenes$Gene)
    }
    
    
    if(length(intersect(cid,colnames(RNAseq.fpkms)))>0){
        
        expProf<-RNAseq.fpkms[,cid]
        notExpGenes<-names(which(expProf<th))
        
        if(amplified){
            notExpGenes<-intersect(notExpGenes,amplifiedGenes)
        }
        return(notExpGenes)
    }else{
        print('No RNAseq data available for this cell line')
        return(NULL)
    }
}
ccr.get.CCLEgisticSets<-function(cellLine,CCLE.gisticCNA=NULL){
    
    data(GDSC.CL_annotation,envir = environment())
    
    if(is.null(CCLE.gisticCNA)){
        data(CCLE.gisticCNA,envir = environment())    
    }
    
    
    if(!is.element(cellLine,GDSC.CL_annotation$CL.name) & !is.element(cellLine,GDSC.CL_annotation$COSMIC.ID)){
        stop('Cell line not found: Please provide a valid cell line cosmic identifier or name (see available cell lines here: http://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources//Data/suppData/TableS1E.xlsx)', call. = TRUE, domain = NULL)
    }
    
    if(is.element(cellLine,GDSC.CL_annotation$COSMIC.ID)){
        cid<-cellLine
    }else{
        cid<-as.character(GDSC.CL_annotation$COSMIC.ID[GDSC.CL_annotation$CL.name==cellLine])
    }
    
    if(length(intersect(cid,colnames(CCLE.gisticCNA)))>0){
        gisticProf<-CCLE.gisticCNA[,cid]
        
        tmpGe<-rownames(CCLE.gisticCNA)[which(gisticProf== -2)]
        gistic_min2<-tmpGe
        
        tmpGe<-rownames(CCLE.gisticCNA)[which(gisticProf== -1)]
        gistic_min1<-tmpGe
        
        tmpGe<-rownames(CCLE.gisticCNA)[which(gisticProf== 0)]
        gistic_wt<-tmpGe
        
        tmpGe<-rownames(CCLE.gisticCNA)[which(gisticProf== 1)]
        gistic_plus1<-tmpGe
        
        tmpGe<-rownames(CCLE.gisticCNA)[which(gisticProf== 2)]
        gistic_plus2<-tmpGe
        
        return(list(gm2=gistic_min2,gm1=gistic_min1,
                    gz=gistic_wt,
                    gp1=gistic_plus1,gp2=gistic_plus2))
    }else{
        print('No gistic CNA scores available for this cell line')
        return(NULL)
    }
    
}


#### Assessment and visualisation
ccr.PrecisionRecallCurve<-function(FCsprofile,positives,negatives,display=TRUE,FDRth=NULL){
    
    FCsprofile<-FCsprofile[intersect(c(positives,negatives),names(FCsprofile))]
    
    predictions<-FCsprofile
    observations<-is.element(names(FCsprofile),positives)+0
    names(observations)<-names(predictions)
    
    RES<-roc(observations,predictions,direction = '>')
    
    if(display){
        plot(RES,col='blue',lwd=3,xlab='TNR',ylab='Recall')
    }
    
    SENS<-NULL
    threshold<-NULL
    if(length(FDRth)>0){
        COORS<-coords(RES,'all',ret = c('threshold','ppv','sensitivity','specificity'))
        
        FDR5percTh<-max(COORS['threshold',which(COORS['ppv',]>=(1-FDRth))])
        
        threshold<-COORS['threshold',min(which(COORS['threshold',]<=FDR5percTh))]
        
        SENS<-COORS['sensitivity',min(which(COORS['threshold',]<=FDR5percTh))]
        SPEC<-COORS['specificity',min(which(COORS['threshold',]<=FDR5percTh))]
        if(display){
            abline(h=SENS,lty=2)
        }
    }
    
    if(display){
        if(length(SENS)==0){
            legend('bottomright',paste('AUC = ',format(RES$auc,digits=3)),bty = 'n')    
        }else{
            legend('bottomright',c(paste('Recall ',100*FDRth,'%FDR = ',format(SENS,digits=3),sep=''),
                                   paste('AUC = ',format(RES$auc,digits=3))),bty = 'n')
        }
        
    }    
    
    RES<-list(AUC=RES$auc,Recall=SENS,sigthreshold=threshold)
    ### threshold, and recall at fixed FDR to be returned
    return(RES)
}
ccr.VisDepAndSig<-function(FCsprofile,
                           SIGNATURES,
                           TITLE='',
                           pIs=NULL,
                           nIs=NULL,
                           th=0.05,
                           plotFCprofile=TRUE){
    
    sigNames<-names(SIGNATURES)
    par(mar=c(5,5,8,0))
    
    nsig<-length(SIGNATURES)
    
    
    if(length(pIs)>0 & length(nIs)>0){
        presentGenes<-intersect(c(SIGNATURES[[pIs]],SIGNATURES[[nIs]]),names(FCsprofile))
        predictions<-FCsprofile[presentGenes]
        observations<-is.element(presentGenes,SIGNATURES[[pIs]])+0
        names(observations)<-presentGenes
        RES<-roc(observations,predictions,'>')
        COORS<-coords(RES,'all',ret = c('threshold','ppv'))
        FDR5percTh<-max(COORS['threshold',which(COORS['ppv',]>=(1-th))])
        FDR5percRANK<-max(which(sort(FCsprofile)<=FDR5percTh))
    }else{
        FDR5percRANK<-NULL
    }
    
    layout(t(matrix(c(rep(1,5),1:(nsig+1),rep(1,5),1:(nsig+1)),nsig+6,2)))
    
    nelements<-length(FCsprofile)
    
    if (plotFCprofile){
        
        plot(sort(FCsprofile,decreasing=TRUE),length(FCsprofile):1,
             ylim=c(nelements,1),pch=16,frame.plot = FALSE,yaxt='n',xlim=c(min(FCsprofile),max(FCsprofile)+1),
             ylab='Depletion Rank',xlab='Log FC',log='y',cex=1.2,cex.lab=1.5,cex.axis=1.2,main=TITLE,cex.main=1.5)
    }else{
        plot(0,0,
             ylim=c(nelements,1),pch=16,frame.plot = FALSE,yaxt='n',xlim=c(min(FCsprofile),max(FCsprofile)+1),
             ylab='Depletion Rank',xlab='Log FC',log='y',cex=1.2,cex.lab=1.5,cex.axis=1.2,main=TITLE,cex.main=1.5)
    }
    lines(x=c(0,0),y=c(1,(length(FCsprofile))+10000),lty=2)
    
    axis(side=2,c(1,10,100,1000,10000,nelements),
         labels = c(1,10,100,1000,10000,nelements),
         las=2,cex.axis=1)
    par(xpd=TRUE)
    lines(c(min(FCsprofile),max(FCsprofile)+1),c(FDR5percRANK,FDR5percRANK),
          col='red',lwd=3,lty=2)
    
    text(x=max(FCsprofile),y=FDR5percRANK,'5%FDR',pos = 3,col='red',cex=1)
    
    TPR<-vector()
    
    N<-length(FCsprofile)
    n<-FDR5percRANK
    
    FTP<-vector()
    PPV<-vector()
    par(xpd=TRUE)
    for (i in 1:length(SIGNATURES)){
        par(mar=c(5,0,8,0))
        
        plot(1,1,xlim=c(0,1),
             ylim=c(length(FCsprofile),1),
             xaxt='n',yaxt='n',ylab='',xlab='',col=NA,log='y',frame.plot=FALSE)
        
        hitPositions<-match(SIGNATURES[[i]],names(sort(FCsprofile)))
        hitPositions<-hitPositions[!is.na(hitPositions)]
        
        abline(h=hitPositions[hitPositions<=FDR5percRANK],lwd=2,col='blue')
        abline(h=hitPositions[hitPositions>FDR5percRANK],lwd=2,col='gray')
        
        
        TPR[i]<-length(which(hitPositions<=FDR5percRANK))/length(hitPositions)
        
        lines(c(0,1),c(FDR5percRANK,FDR5percRANK),col='red',lwd=3,lty=2)
        
        text(0.4,1,labels = sigNames[i],pos = 4,offset = 0,srt=80,cex=1)
        
        if(i<(length(SIGNATURES)-1)){
            text(0.6,nelements+50000,paste(round(100*TPR[i]),'%',sep=''),
                 cex=1.2,col='blue')
        }
        
    }
    
    names(TPR)<-sigNames
    
    return(TPR)
    
}

#### Utils not documented
#### Assessment and visualisation not documented



ccr.assessmentReport<-function(cellLine,
                                libraryAnnotation,
                                correctedFCs,
                                sgRNAlevel=TRUE,
                                outDir='.',
                               GDSC.geneLevCNA=NULL,
                               CCLE.gisticCNA=NULL,
                               RNAseq.fpkms=NULL){
    
    
    if(is.null(GDSC.geneLevCNA)){
        data(GDSC.geneLevCNA)
    }
    if(is.null(CCLE.gisticCNA)){
        data(CCLE.gisticCNA)
    }
    if(is.null(RNAseq.fpkms)){
        data(RNAseq.fpkms)
    }
    
    
    GDSC.cna<-ccr.get.gdsc1000.AMPgenes(cellLine = cellLine,GDSC.geneLevCNA = GDSC.geneLevCNA,minCN = 0)
    CCLE.cna<-ccr.get.CCLEgisticSets(cellLine = cellLine,CCLE.gisticCNA = CCLE.gisticCNA)
    notExp<-ccr.get.nonExpGenes(cellLine = cellLine,RNAseq.fpkms = RNAseq.fpkms)
    
    data(BAGEL_essential)
    data(BAGEL_nonEssential)
    
    
    
    nnames<-c('Dep (PN)',
              'Dep (Gistic)',
              'notExp',
              'Amp (Gistic +1)',
              'Amp (Gistic +2)',
              'Amp (Gistic +1) notExp',
              'Amp (Gistic +2) notExp',
              'Amp (PNs >= 2)',
              'Amp (PNs >= 4)',
              'Amp (PNs >= 8)',
              'Amp (PNs >= 10)',
              'Amp (PNs >= 2) notExp',
              'Amp (PNs >= 4) notExp',
              'Amp (PNs >= 8) notExp',
              'Amp (PNs >= 10) notExp',
              'Essential',
              'BAGEL Essential',
              'BAGEL Ess.Only',
              'BAGEL nonEssential',
              'whole-library')
    
    
}
# 
#     whatToPlot<-c('avgFC','correctedFC')
    
    # #pdf(paste(outDir,cellLine,'_CVgistic_densities.pdf',sep=''),width = 5,height = 8.5)
    # 
    # par(mfrow=c(2,1))
    # 
    # Xmin<-min(min(correctedFCs[,'avgFC']),
    #           min(correctedFCs[,'correctedFC']))-1
    # Xmax<-max(max(correctedFCs[,'avgFC']),
    #           max(correctedFCs[,'correctedFC']))+1
    # 
    # for (j in 1:2){
    #     
    #     if(length(gisticGuides$gp1)>0){
    #         densGistAMP1<-
    #             density(correctedFCs[gisticSet$gp1,whatToPlot[j]],na.rm = TRUE)
    #     }else{
    #         densGistAMP1<-list(x=NA,y=NA)
    #     }
    #     
    #     if(length(gisticGuides$gp2)>0){
    #         densGistAMP2<-density(correctedFCs$corrected_logFCs[gisticGuides$gp2,whatToPlot[j]],na.rm = TRUE)
    #     }else{
    #         densGistAMP2<-list(x=NA,y=NA)
    #     }
    #     
    #     densOthers<-density(correctedFCs$corrected_logFCs[setdiff(1:nrow(correctedFCs$corrected_logFCs),
    #                                                               c(gisticGuides$gp1,gisticGuides$gp2)),
    #                                                       whatToPlot[j]],na.rm = TRUE)
    #     
    #     COLS<-c('red2','red4','darkgrey')
    #     toPlot<-list(densGistAMP1,densGistAMP2,densOthers)
    #     
    #     if (j==1){
    #         xlab=''
    #         title='pre CRISPRcleanR'
    #     }else{
    #         xlab='sgRNA log FC'
    #         title='post CRISPRcleanR'
    #     }
    #     if(i==1){
    #         par(mar=c(2,4,4,1))
    #     }else{
    #         par(mar=c(4,4,2,1))
    #     }
    #     
    #     ccr.multDensPlot(TOPLOT = toPlot,COLS = COLS,XLIMS = c(Xmin,Xmax),TITLE=title,
    #                      LEGentries = c('Amp (Gistic +1)','Amp (Gistic +2)','others'),XLAB = xlab)
    # }
    # 
    # dev.off()
    
    
    # pdf(paste(outDir,cellLine,'_CVpicnic_densities.pdf',sep=''),width = 5,height = 8.5)
    # 
    # par(mfrow=c(2,1))
    # 
    # for (j in 1:2){
    #     
    #     
    #     currentG<-ampGenes$Gene[which(ampGenes$CN>=8)]
    #     
    #     if(sgRNAlevel){
    #         currentG<-ccr.genes2sgRNAs(libraryAnnotation,currentG)    
    #     }
    #     
    #     if(length(currentG)>0){
    #         densPicNicAMP8<-
    #             density(correctedFCs[currentGuides,whatToPlot[j]],na.rm = TRUE)
    #     }else{
    #         densPicNicAMP8<-list(x=NA,y=NA)
    #     }
    #     
    #     currentG<-ampGenes$Gene[which(ampGenes$CN>=10)]
    #     if(sgRNAlevel){
    #         currentG<-ccr.genes2sgRNAs(libraryAnnotation,currentG)    
    #     }
    #     if(length(currentG)>0){
    #         densPicNicAMP10<-density(correctedFCs[currentG,whatToPlot[j]],na.rm = TRUE)
    #     }else{
    #         densPicNicAMP10<-list(x=NA,y=NA)
    #     }
    #     
    #     
    #     densOthers<-density(correctedFCs[setdiff(1:nrow(correctedFCs),
    #                                                               c(densPicNicAMP10,densPicNicAMP8)),
    #                                                       whatToPlot[j]],na.rm = TRUE)
    #     
    #     COLS<-c("#FF243C","#FF0000",'darkgrey')
    #     toPlot<-list(densPicNicAMP8,densPicNicAMP10,densOthers)
    #     
    #     if (j==1){
    #         xlab=''
    #         title='pre CRISPRcleanR'
    #     }else{
    #         xlab='sgRNA log FC'
    #         title='post CRISPRcleanR'
    #     }
    #     if(i==1){
    #         par(mar=c(2,4,4,1))
    #     }else{
    #         par(mar=c(4,4,2,1))
    #     }
    #     
    #     ccr.multDensPlot(TOPLOT = toPlot,
    #                      COLS = COLS,XLIMS = c(Xmin,Xmax),TITLE=title,
    #                      LEGentries = c('Amp (PNs >= 8)','Amp (PNs >= 10)','Others'),XLAB = xlab)
    # }
    # 
    # dev.off()
    # 
    # pdf(paste(outDir,cellLine,'_CVE_densities.pdf',sep=''),width = 5,height = 8.5)
    # par(mfrow=c(2,1))
    # 
    # for (j in 1:2){
    #     
    #     notExpAmpGuidesp8<-intersect(names(which(ampGuidesPattern$MINcna>=8)),notExpGuides)
    #     
    #     if(length(notExpAmpGuidesp8)>0){
    #         densNotExpAMPp8<-density(correctedFCs$corrected_logFCs[notExpAmpGuidesp8,whatToPlot[j]],na.rm = TRUE)
    #     }else{
    #         densNotExpAMPp8<-list(x=NA,y=NA)
    #     }
    #     
    #     notExpAmpGuidesp10<-intersect(names(which(ampGuidesPattern$MINcna>=10)),notExpGuides)
    #     
    #     if(length(notExpAmpGuidesp10)>0){
    #         densNotExpAMPp10<-density(correctedFCs$corrected_logFCs[notExpAmpGuidesp10,whatToPlot[j]],na.rm = TRUE)
    #     }else{
    #         densNotExpAMPp10<-list(x=NA,y=NA)
    #     }
    #     
    #     notExpAmpG<-intersect(gisticGuides$gp2,notExpGuides)
    #     
    #     if(length(notExpAmpG)>0){
    #         densNotExpAmpG<-density(correctedFCs$corrected_logFCs[notExpAmpG,whatToPlot[j]],na.rm = TRUE)
    #     }else{
    #         densNotExpAmpG<-list(x=NA,y=NA)
    #     }
    #     
    #     densOthers<-density(correctedFCs$corrected_logFCs[setdiff(1:nrow(correctedFCs$corrected_logFCs),
    #                                                               c(notExpAmpGuidesp8,notExpAmpGuidesp10,notExpAmpG)),
    #                                                       whatToPlot[j]],na.rm = TRUE)
    #     
    #     COLS<-c('blue',"#2F4CAC","#00008B",'darkgray')
    #     toPlot<-list(densNotExpAmpG,densNotExpAMPp8,densNotExpAMPp10,densOthers)
    #     
    #     if (j==1){
    #         xlab=''
    #         title='pre CRISPRcleanR'
    #     }else{
    #         xlab='sgRNA log FC'
    #         title='post CRISPRcleanR'
    #     }
    #     if(i==1){
    #         par(mar=c(2,4,4,1))
    #     }else{
    #         par(mar=c(4,4,2,1))
    #     }
    #     
    #     ccr.multDensPlot(TOPLOT = toPlot,COLS = COLS,XLIMS = c(Xmin,Xmax),TITLE=title,
    #                      LEGentries = c('Gistic +2 notExp','PNs >= 8 notExp','PNs >= 10 notExp','Others'),XLAB = xlab)
    #     
    # }
    # 
    # dev.off()
    # 
    # 
    # pdf(paste(outDir,cellLine,'_ENE_densities.pdf',sep=''),width = 5,height = 8.5)
    # par(mfrow=c(2,1))
    # for (j in 1:2){
    #     densEssGenes<-density(correctedFCs$corrected_logFCs[essGeneGuides,whatToPlot[j]],na.rm = TRUE)
    #     densBagelEssGenes<-density(correctedFCs$corrected_logFCs[essBagelGeneGuides,whatToPlot[j]],na.rm = TRUE)
    #     densBagelOnlyEssGenes<-density(correctedFCs$corrected_logFCs[setdiff(essBagelGeneGuides,essGeneGuides),whatToPlot[j]],na.rm = TRUE)
    #     densBagelNonEssGenes<-density(correctedFCs$corrected_logFCs[NONessBagelGeneGuides,whatToPlot[j]],na.rm=TRUE)
    #     densNonExpGenes<-density(correctedFCs$corrected_logFCs[NONessBagelGeneGuides,whatToPlot[j]],na.rm=TRUE)
    #     
    #     COLS<-c('darkgreen','green','darkcyan',
    #             "bisque4")
    #     toPlot<-list(densEssGenes,densBagelEssGenes,densBagelOnlyEssGenes,densBagelNonEssGenes)
    #     
    #     if (j==1){
    #         xlab=''
    #         title='pre CRISPRcleanR'
    #     }else{
    #         xlab='sgRNA log FC'
    #         title='post CRISPRcleanR'
    #     }
    #     if(i==1){
    #         par(mar=c(2,4,4,1))
    #     }else{
    #         par(mar=c(4,4,2,1))
    #     }
    #     
    #     ccr.multDensPlot(TOPLOT = toPlot,COLS = COLS,XLIMS = c(Xmin,Xmax),TITLE=title,
    #                      LEGentries = c('Essential','BAGEL Essential','BAGEL Ess.Only','BAGEL nonEssential'),XLAB = xlab)
    # }
    # 
    # dev.off()
    # 
    # CS<-list(essGeneGuides,essBagelGeneGuides,essBagelOnlyGuides,
    #          names(which(ampGuidesPattern$MINcna>=8)),
    #          gisticGuides$gp2,
    #          notExpGuides,
    #          notExpAmpGuidesp8,
    #          NONessBagelGeneGuides)
    # 
    # pdf(paste(outDir,cellLine,'_ROCs.pdf',sep=''),width = 8,height = 14)
    # 
    # COLS<-c('darkgreen','green','darkcyan',
    #         "red",
    #         "red4",
    #         "#626262",
    #         'blue',
    #         "bisque4")
    # 
    # par(mfrow=c(2,1))
    # par(mar=c(4,4,4,1))
    # O1<-order(correctedFCs$corrected_logFCs$avgFC)
    # predictions1<-rownames(correctedFCs$corrected_logFCs)[O1]
    # 
    # 
    # AUCs<-vector()
    # 
    # dd<-c(2,4,5,6,7,8)
    # for (i in 1:length(CS)){
    #     
    #     currentSet<-CS[[i]]
    #     currentSet<-intersect(currentSet,predictions1)
    #     pp1<-cumsum(is.element(predictions1,currentSet))/length(currentSet)
    #     
    #     if(i == 1){
    #         plot(100*(1:length(predictions1))/length(predictions1),100*pp1,type='l',xlim=c(0,100),ylim=c(0,100),ylab='cumulative %',
    #              col=COLS[i],xlab='',
    #              main=cellLine,lwd=2)
    #     }else{
    #         lines(100*(1:length(predictions1))/length(predictions1),100*pp1,type='l',xlim=c(0,100),ylim=c(0,100),col=COLS[i],lwd=2)
    #     }
    #     
    #     AUCs[i]<-trapz(1:length(predictions1)/length(predictions1),pp1)
    # }
    # 
    # legend('bottomright',legend = c('Essential',
    #                                 'Bagel Ess.',
    #                                 'Bagel Ess.only',
    #                                 'Amp.Genes (PN >= 8)',
    #                                 'Amp.Genes (Gistic = 2)',
    #                                 'NotExpAmplified',
    #                                 'NotExp',
    #                                 'Non.Ess'),col=COLS[c(1,2,3,4,5,7,6,8)],lty=1,cex = 1,lwd=2)
    # 
    # O2<-order(correctedFCs$corrected_logFCs$correctedFC)
    # predictions2<-rownames(correctedFCs$corrected_logFCs)[O2]
    # 
    # AUCs1<-vector()
    # for (i in 1:length(CS)){
    #     currentSet<-CS[[i]]
    #     currentSet<-intersect(currentSet,predictions2)
    #     
    #     pp2<-cumsum(is.element(predictions2,currentSet))/length(currentSet)
    #     
    #     if(i == 1){
    #         plot(100*(1:length(predictions2))/length(predictions2),100*pp2,type='l',xlim=c(0,100),ylim=c(0,100),ylab='cumulative %',col=COLS[i],
    #              xlab='sgRNA FC percentile',main=paste(cellLine,'post CRISPRcleanR'),lwd=2)
    #     }else{
    #         lines(100*(1:length(predictions2))/length(predictions2),100*pp2,type='l',xlim=c(0,100),ylim=c(0,100),col=COLS[i],lwd=2)
    #     }
    #     
    #     AUCs1[i]<-trapz(1:length(predictions2)/length(predictions2),pp2)
    # }
    # 
    # 
    # legend('bottomright',legend = c('Essential',
    #                                 'Bagel Ess.',
    #                                 'Bagel Ess.only',
    #                                 'Amp.Genes (PN >= 8)',
    #                                 'Amp.Genes (Gistic = 2)',
    #                                 'NotExpAmplified',
    #                                 'NotExp',
    #                                 'Non.Ess'),col=COLS[c(1,2,3,4,5,7,6,8)],lty=1,cex = 1,lwd=2)
    # dev.off()
    # 
    # AUCs<-rbind(AUCs,AUCs1)
    # rownames(AUCs)<-c('pre-CRISPRcleanR','post-CRISPRcleanR')
    # colnames(AUCs)<-c('Ess.Genes','Bagel.Ess','Bagel.Ess.Only','Amp.Genes (PN>=8)','Amp.Genes (Gistic=2)','NotExp','NotExpAmp','Non.Ess')
    # 
    # EFFsizes<-rbind(EFFsizes,EFFsizesN)
    # rownames(EFFsizes)<-c('pre-CRISPRcleanR','post-CRISPRcleanR')
    # colnames(EFFsizes)<-nnames[1:19]
    # 
    # PVALS<-rbind(PVALS,PVALSn)
    # rownames(EFFsizes)<-c('pre-CRISPRcleanR','post-CRISPRcleanR')
    # colnames(EFFsizes)<-nnames[1:19]
    # 
    # SIGNS<-rbind(SIGNS,SIGNSn)
    # rownames(SIGNS)<-c('pre-CRISPRcleanR','post-CRISPRcleanR')
    # colnames(SIGNS)<-nnames[1:19]
    # 
    # 
    # idx<-c(16,17,18,10,7,3,14,19)
    # EFFsizes<-EFFsizes[,idx]
    # colnames(EFFsizes)<-colnames(AUCs)
    # 
    # PVALS<-PVALS[,idx]
    # colnames(PVALS)<-colnames(AUCs)
    # 
    # SIGNS<-SIGNS[,idx]
    # colnames(SIGNS)<-colnames(AUCs)
    # 
    # return(list(AUCs=AUCs,PVALS=PVALS,SIGNS=SIGNS,EFFsizes=EFFsizes,preCorrValues=preValues,postCorrValues=postValues))
#}


## other exported functions

## not exported functions
ccr.boxplot<-function(toPlot,main,names){
         
         boxplot(toPlot,main=main,names = names,las=2)
         MEANS<-apply(toPlot,MARGIN = 2,FUN = 'mean')
         SD<-apply(toPlot,MARGIN = 2,FUN = 'sd')
         for (i in 1:length(MEANS)){
             lines(x=c(i-0.3,i+0.3),y=c(MEANS[i],MEANS[i])+SD[i],col='blue',lwd=2)
             lines(x=c(i-0.2,i+0.2),y=c(MEANS[i],MEANS[i]),col='red',lwd=5)
             lines(x=c(i-0.3,i+0.3),y=c(MEANS[i],MEANS[i])-SD[i],col='blue',lwd=2)
         }
    }

