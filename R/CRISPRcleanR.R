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
    
    set.seed(0xA5EED)
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
    
    
    
    adjusted<-
        as.data.frame(cbind(normalised_counts$sgRNA,
                            normalised_counts$gene,
                            correctedCounts),stringAsFactors=FALSE)
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
    
    data(GDSC.CL_annotation,envir = environment()) 
    
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
ccr.PlainTsvFile<-function(sgRNA_count_object,fprefix='',
                           path='./'){
    
    currentFileContent<-sgRNA_count_object
    
    fname<-paste(path,fprefix,'_sgRNA_count.tsv',sep='')
    
    write.table(currentFileContent,quote=FALSE,row.names = FALSE,sep='\t',
                file=fname)
    
    return(fname)
}

#### Assessment and visualisation
ccr.ROC_Curve<-function(FCsprofile,positives,negatives,display=TRUE,FDRth=NULL){
    
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
    COORS<-coords(RES,'all',ret = c('threshold','ppv','sensitivity','specificity'))
    if(length(FDRth)>0){
        
        
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
    
    
    COORS<-t(COORS[c('specificity','sensitivity','threshold'),])
    RES<-list(AUC=RES$auc,Recall=SENS,sigthreshold=threshold,curve=COORS)
    ### threshold, and recall at fixed FDR to be returned
    return(RES)
}
ccr.PrRc_Curve<-function(FCsprofile,positives,negatives,display=TRUE,FDRth=NULL){
    
    FCsprofile<-FCsprofile[intersect(c(positives,negatives),names(FCsprofile))]
    
    predictions<- -FCsprofile
    observations<-is.element(names(FCsprofile),positives)+0
    names(observations)<-names(predictions)
    
    
    prc<-pr.curve(scores.class0 = predictions,weights.class0 = observations,
                  curve = TRUE,sorted = TRUE)

    PRECISION<-prc$curve[,2]
    RECALL<-prc$curve[,1]
    
    if(display){
        plot(RECALL,PRECISION,col='blue',lwd=3,xlab='Recall',ylab='Precision',type='l',xlim=c(0,1),ylim=c(0,1))
    }
    
    SENS<-NULL
    threshold<-NULL
    if(length(FDRth)>0){
        
        FDR5percTh<- -prc$curve[min(which(prc$curve[,2]>= 1-FDRth)),3]
        SENS<- prc$curve[min(which(prc$curve[,2]>= 1-FDRth)),1]
        threshold<-FDR5percTh
        if(display){
            abline(h=1-FDRth,lty=2)
            
            abline(v=SENS,lty=1)
        }
    }
    
    if(display){
        if(length(SENS)==0){
            legend('bottomleft',paste('AUC = ',format(prc$auc.integral,digits=3)),bty = 'n')    
        }else{
            legend('bottomleft',c(paste('Recall ',100*FDRth,'%FDR = ',format(SENS,digits=3),sep=''),
                                  paste('AUC = ',format(prc$auc.integral,digits=3))),bty = 'n')
        }
        
        abline(h=sum(observations)/length(observations))
    }    
    # 
    curve<-prc$curve
    colnames(curve)<-c('recall','precision','threshold')
    RES<-list(AUC=prc$auc.integral,Recall=SENS,sigthreshold=threshold,curve=curve)
    # ### threshold, and recall at fixed FDR to be returned
    return(RES)
}
ccr.randomised_ROC<-function(FCs,PERCrandn,ntrials,positives,negatives,LibraryAnnotation){
    
    
    posGuides<-ccr.genes2sgRNAs(libraryAnnotation = LibraryAnnotation,genes = positives)
    negGuides<-ccr.genes2sgRNAs(libraryAnnotation = LibraryAnnotation,genes = negatives)
    
    guidesToShuffle<-intersect(union(posGuides,negGuides),names(FCs))
    
    idGuidesToShuffle<-match(guidesToShuffle,names(FCs))
    
    nguides<-round(length(guidesToShuffle)*PERCrandn/100)
    
    rndSENSITIVITY <- NULL
    rndSPECIFICITY <- NULL
    rndAUROC <- NULL
    
    rndRECALL <- NULL
    rndPRECISION <- NULL
    rndAUPRC <- NULL
    
    rndRECALL_at_fixedFDR<-NULL
    
    
    
    for (i in 1:ntrials){
        print(i)
        
        RNDFCs<-FCs
        mid<-sample(length(guidesToShuffle),nguides)
        toShuffle<-idGuidesToShuffle[mid]    
        
        otherGuides<-setdiff(1:length(FCs),toShuffle)
        
        toReplace <- otherGuides[sample(length(otherGuides), nguides)]
        
        names(RNDFCs)[c(toShuffle,toReplace)] <- names(FCs)[c(toReplace,toShuffle)]
        
        RNDgeneFCs <- ccr.geneMeanFCs(RNDFCs, LibraryAnnotation)
        
        RESroc <-
            ccr.ROC_Curve(
                RNDgeneFCs,
                positives,
                negatives,
                FDRth = 0.05,
                display = FALSE
            )
        rndSENSITIVITY <- rbind(rndSENSITIVITY, RESroc$curve[, 'sensitivity'])
        rndSPECIFICITY <- rbind(rndSPECIFICITY, RESroc$curve[, 'specificity'])
        rndAUROC <- c(rndAUROC, RESroc$AUC)
        
        RESprrc <-
            ccr.PrRc_Curve(
                RNDgeneFCs,
                positives,
                negatives,
                FDRth = 0.05,
                display = FALSE
            )
        rndRECALL <- rbind(rndRECALL, RESprrc$curve[, 'recall'])
        rndPRECISION <- rbind(rndPRECISION, RESprrc$curve[, 'precision'])
        rndAUPRC <- c(rndAUPRC, RESprrc$AUC)
        rndRECALL_at_fixedFDR<-c(rndRECALL_at_fixedFDR,RESprrc$Recall)
    }
    
    return(list(rndSENSITIVITY = rndSENSITIVITY,
                rndSPECIFICITY = rndSPECIFICITY,
                rndAUROC = rndAUROC,
                rndRECALL = rndRECALL,
                rndPRECISION = rndPRECISION,
                rndAUPRC = rndAUPRC,
                rndRECALL_at_fixedFDR = rndRECALL_at_fixedFDR))
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
        RES<-ccr.fixedFDRthreshold(FCsprofile = FCsprofile,
                                   TruePositives = SIGNATURES[[pIs]],
                                   TrueNegatives = SIGNATURES[[nIs]],
                                   th = th)
        FDR5percRANK<-RES$RANK
        FDR5percTh<-RES$FCth
        
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

ccr.perf_statTests<-function(cellLine,
                                libraryAnnotation,
                                correctedFCs,
                                outDir='./',
                               GDSC.geneLevCNA=NULL,
                               CCLE.gisticCNA=NULL,
                               RNAseq.fpkms=NULL){
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
        
        guideSets<-ccr.get.guideSets(cellLine,GDSC.geneLevCNA,CCLE.gisticCNA,RNAseq.fpkms,libraryAnnotation = libraryAnnotation)
        PVALS<-vector()
        PVALSn<-vector()
        
        EFFsizes<-vector()
        EFFsizesN<-vector()
        
        SIGNS<-vector()
        SIGNSn<-vector()
    
        for (i in 1:length(guideSets)){
            
            print(paste('Testing sgRNAs targeting:', nnames[i], 'genes'))
            currentSet<-guideSets[[i]]
            
            if(length(currentSet)>1){
                tt<-t.test(correctedFCs[currentSet,'avgFC'],
                           correctedFCs[setdiff(rownames(correctedFCs),currentSet),'avgFC'])
                
                ttn<-t.test(correctedFCs[currentSet,'correctedFC'],
                            correctedFCs[setdiff(rownames(correctedFCs),currentSet),'correctedFC'])
                
                EFFsizes[i]<-ccr.cohens_d(x = correctedFCs[currentSet,'avgFC'],
                                          y = correctedFCs[setdiff(rownames(correctedFCs),currentSet),'avgFC'])
                EFFsizesN[i]<-ccr.cohens_d(x = correctedFCs[currentSet,'correctedFC'],
                                           y = correctedFCs[setdiff(rownames(correctedFCs),currentSet),'correctedFC'])
                
                PVALS[i]<-tt$p.value
                PVALSn[i]<-ttn$p.value
                
                SIGNS[i]<-sign(tt$estimate[2]-tt$estimate[1])
                SIGNSn[i]<-sign(ttn$estimate[2]-ttn$estimate[1])
                
            }else{
                PVALS[i]<-NA
                PVALSn[i]<-NA
                
                EFFsizes[i]<-NA
                EFFsizesN[i]<-NA
                
                SIGNS[i]<-NA
                SIGNSn[i]<-NA
                
            }
        }
        
        pdf(paste(outDir,cellLine,'_bp.pdf',sep=''),width = 15,height = 10)
        AT<-c(1,2,3,4.5,5.5,6.5,7.5,9,10,11,12,13,14,15,16,17.5,18.5,19.5,20.5,22)
        par(mfrow=c(2,1))
        par(mar=c(6,4,4,1))
        par(xpd=NA)
        plotpar<-boxplot(correctedFCs[guideSets[[1]],'avgFC'],
                         correctedFCs[guideSets[[2]],'avgFC'],
                         correctedFCs[guideSets[[3]],'avgFC'],
                         correctedFCs[guideSets[[4]],'avgFC'],
                         correctedFCs[guideSets[[5]],'avgFC'],
                         correctedFCs[guideSets[[6]],'avgFC'],
                         correctedFCs[guideSets[[7]],'avgFC'],
                         correctedFCs[guideSets[[8]],'avgFC'],
                         correctedFCs[guideSets[[9]],'avgFC'],
                         correctedFCs[guideSets[[10]],'avgFC'],
                         correctedFCs[guideSets[[11]],'avgFC'],
                         correctedFCs[guideSets[[12]],'avgFC'],
                         correctedFCs[guideSets[[13]],'avgFC'],
                         correctedFCs[guideSets[[14]],'avgFC'],
                         correctedFCs[guideSets[[15]],'avgFC'],
                         correctedFCs[guideSets[[16]],'avgFC'],
                         correctedFCs[guideSets[[17]],'avgFC'],
                         correctedFCs[guideSets[[18]],'avgFC'],
                         correctedFCs[guideSets[[19]],'avgFC'],
                         correctedFCs$avgFC,
                         at=AT,names = nnames,las=2,outline = FALSE,
                         frame.plot=FALSE,main=cellLine,ylab='pre-CRISPRcleanR sgRNA log2 FC',
                         col=c('#B0B0B0',"#898989","#626262",c('red2','red4'),'cornflowerblue','blue',
                               c("#FF6EB4","#FF4978","#FF243C","#FF0000"),"#8EE5EE","#5E98CD","#2F4CAC","#00008B",'darkgreen','green','darkcyan',
                               "bisque4","white"))
        
        sigVec<-rep('',19)
        sigVec[which(PVALS<0.05)]<-'.'
        sigVec[which(PVALS<0.05 & EFFsizes > 0.5)]<-'*'
        sigVec[which(PVALS<0.05 & EFFsizes > 1)]<-'**'
        sigVec[which(PVALS<0.05 & EFFsizes > 2)]<-'***'
        
        posY<-(max(correctedFCs$avgFC)-min(correctedFCs$avgFC))/30+max(correctedFCs$avgFC)
        
        for (i in 1:19){
            text(AT[i],posY,sigVec[i],cex = 1.5)
        }
        dd<-mean(correctedFCs$avgFC)
        lines(x = c(0.2,22.8),y = c(dd,dd),lty=2,col='red',lwd=2)
        
        #text(AT+0.25,min(correctedFCs$corrected_logFCs$avgFC)-1,labels = nnames,srt=45,pos = 2)
        par(mar=c(2,4,8,1))
        par(xpd=NA)
        
        boxplot(correctedFCs[guideSets[[1]],'correctedFC'],
                correctedFCs[guideSets[[2]],'correctedFC'],
                correctedFCs[guideSets[[3]],'correctedFC'],
                correctedFCs[guideSets[[4]],'correctedFC'],
                correctedFCs[guideSets[[5]],'correctedFC'],
                correctedFCs[guideSets[[6]],'correctedFC'],
                correctedFCs[guideSets[[7]],'correctedFC'],
                correctedFCs[guideSets[[8]],'correctedFC'],
                correctedFCs[guideSets[[9]],'correctedFC'],
                correctedFCs[guideSets[[10]],'correctedFC'],
                correctedFCs[guideSets[[11]],'correctedFC'],
                correctedFCs[guideSets[[12]],'correctedFC'],
                correctedFCs[guideSets[[13]],'correctedFC'],
                correctedFCs[guideSets[[14]],'correctedFC'],
                correctedFCs[guideSets[[15]],'correctedFC'],
                correctedFCs[guideSets[[16]],'correctedFC'],
                correctedFCs[guideSets[[17]],'correctedFC'],
                correctedFCs[guideSets[[18]],'correctedFC'],
                correctedFCs[guideSets[[19]],'correctedFC'],
                correctedFCs$correctedFC,
                at=AT,outline = FALSE,
                names = rep('',20),frame.plot=FALSE,main='',ylab='pre-CRISPRcleanR sgRNA log2 FC',
                col=c('#B0B0B0',"#898989","#626262",c('red2','red4'),'cornflowerblue','blue',
                      c("#FF6EB4","#FF4978","#FF243C","#FF0000"),"#8EE5EE","#5E98CD","#2F4CAC","#00008B",'darkgreen','green','darkcyan',
                      "bisque4","white"))
     
        
        sigVec1<-rep('',19)
        sigVec1[which(PVALSn<0.05)]<-'.'
        sigVec1[which(PVALSn<0.05 & EFFsizesN > 0.5)]<-'*'
        sigVec1[which(PVALSn<0.05 & EFFsizesN > 1)]<-'**'
        sigVec1[which(PVALSn<0.05 & EFFsizesN > 2)]<-'***'
        
        posY<-(max(correctedFCs$correctedFC)-min(correctedFCs$correctedFC))/30+max(correctedFCs$correctedFC)
        
        for (i in 1:19){
            text(AT[i],posY,sigVec1[i],cex = 1.5)
        }
        dd<-mean(correctedFCs$correctedFC)
        lines(x = c(0.2,22.8),y = c(dd,dd),lty=2,col='red',lwd=2)
        
        legend('bottomleft',legend = c('.     p < 0.05',
                                       '*     p < 0.05, Effect size > 0.5',
                                       '**    p < 0.05, Effect size > 1',
                                       '***   p < 0.05, Effect size > 2'))
        dev.off()
        
        EFFsizes<-rbind(EFFsizes,EFFsizesN)
        rownames(EFFsizes)<-c('pre-CRISPRcleanR','post-CRISPRcleanR')
        colnames(EFFsizes)<-nnames[1:19]
        
        PVALS<-rbind(PVALS,PVALSn)
        rownames(PVALS)<-c('pre-CRISPRcleanR','post-CRISPRcleanR')
        colnames(PVALS)<-nnames[1:19]
        
        SIGNS<-rbind(SIGNS,SIGNSn)
        rownames(SIGNS)<-c('pre-CRISPRcleanR','post-CRISPRcleanR')
        colnames(SIGNS)<-nnames[1:19]
    
        return(list(PVALS=PVALS,SIGNS=SIGNS,EFFsizes=EFFsizes))
}

ccr.multDensPlot<-function(TOPLOT,COLS,XLIMS,TITLE,LEGentries,XLAB=''){
    
    YM<-vector()
    for (i in 1:length(TOPLOT)){
        YM[i]<-max(TOPLOT[[i]]$y,na.rm = TRUE)
    }
    
    Ymax<-max(YM,na.rm=TRUE)
    
    plot(0,0,col=NA,ylab='density',xlab=XLAB,
         xlim=XLIMS,ylim=c(0,Ymax),type='l',main=TITLE)
    
    for (i in 1:length(TOPLOT)){
        cord.x <- c(TOPLOT[[i]]$x)
        cord.y <- c(TOPLOT[[i]]$y)
        rgbc<-col2rgb(COLS[i])
        currCol<-rgb(rgbc[1],rgbc[2],rgbc[3],alpha = 100,maxColorValue = 255)
        polygon(cord.x,cord.y,col=currCol,border = NA)
        lines(TOPLOT[[i]],col=COLS[i],lwd=3)
        if (i == 1){
            legend('topleft',legend = LEGentries,col=COLS,lwd=3,bty = 'n')
        }
    }
}

ccr.perf_distributions<-function(cellLine,correctedFCs,
                                 GDSC.geneLevCNA=NULL,
                                 CCLE.gisticCNA=NULL,
                                 RNAseq.fpkms=NULL,
                                 minCNs=c(8,10),
                                 libraryAnnotation){
    
    guideSets<-ccr.get.guideSets(cellLine,GDSC.geneLevCNA,CCLE.gisticCNA,RNAseq.fpkms,
                                 libraryAnnotation = libraryAnnotation)
    
    names(guideSets)[16]<-'MSigDB CFEs'
    
    gs<-guideSets
    Xmin<-min(min(correctedFCs[,'avgFC']),min(correctedFCs[,'correctedFC']))-1
    Xmax<-max(max(correctedFCs[,'avgFC']),max(correctedFCs[,'correctedFC']))+1
    
    whatToPlot<-c('avgFC','correctedFC')
    
    geneSetToPlot<-list(AMPLIFIEDpn=c(paste('Amp (PNs >= ',minCNs,')',sep='')),
                        AMPLIFIEDgistic=c("Amp (Gistic +1)","Amp (Gistic +2)"),
                        AMPLIFIEDnotExp=c("Amp (Gistic +2) notExp",paste('Amp (PNs >= ',minCNs,') notExp',sep='')),
                        REFERENCE=c('MSigDB CFEs',"BAGEL Essential","BAGEL nonEssential"))
    
    COLS_list<-list(AMPLIFIEDpn=c("#FF243C","#FF0000",'darkgrey'),
                    AMPLIFIEDgistic=c('red2','red4','darkgrey'),
                    AMPLIFIEDnotExp=c('blue',"#2F4CAC","#00008B",'darkgray'),
                    REFERENCE=c('darkgreen','green',"bisque4"))
    options(warn = -1)
    for (i in 1:4){
        par(mfrow=c(2,1))
        for (j in 1:2){
            
            plotType<-geneSetToPlot[[i]]
            
            densities<-list()
            for (k in 1:length(plotType)){
                if (length(gs[[plotType[k]]]>0)){
                    densities[[k]]<-density(correctedFCs[gs[[plotType[k]]],whatToPlot[j]],na.rm = TRUE)    
                }else{
                    densities[[k]]<-list(x=NA,y=NA)
                }
            }
            
            
            densOthers<-density(correctedFCs[setdiff(rownames(correctedFCs),unlist(gs[plotType])),
                                             whatToPlot[j]],na.rm = TRUE)
            
            COLS<-COLS_list[[i]]
            toPlot<-append(densities,list(densOthers))
            
            if (j==1){
                xlab=''
                title='pre CRISPRcleanR'
            }else{
                xlab='sgRNA log FC'
                title='post CRISPRcleanR'
            }
            
            par(mar=c(4,4,2,1))
            
            options(warn = -1)
            ccr.multDensPlot(TOPLOT = toPlot,COLS = COLS,XLIMS = c(Xmin,Xmax),TITLE=title,
                             LEGentries = c(plotType,'others'),XLAB = xlab)
            
            options(warn = 0)
        }
        options(warn = 0)
    }
    
}

ccr.RecallCurves<-function(cellLine,correctedFCs,
                           GDSC.geneLevCNA=NULL,
                           RNAseq.fpkms=NULL,
                           minCN=8,
                           libraryAnnotation,
                           GeneLev=FALSE){
    
    guideSets<-ccr.get.guideSets(cellLine,GDSC.geneLevCNA,CCLE.gisticCNA=NULL,RNAseq.fpkms,
                                 libraryAnnotation = libraryAnnotation)
    
    if(GeneLev){
        guideSets<-lapply(guideSets,function(x){libraryAnnotation[x,'GENES']})
    }
    
    COLS<-c("bisque4",'green','red','blue')
    
    par(mfrow=c(2,1))
    par(mar=c(4,4,4,1))
    
    toPlot<-c('avgFC','correctedFC')
    
    AUCcombo<-NULL
    for (j in 1:length(toPlot)){
        
        if (j == 1){
            MAIN<-paste(cellLine,'pre-CRISPRcleanR')
        }else{
            MAIN<-paste(cellLine,'post-CRISPRcleanR')
        }
        
        if(!GeneLev){
            O1<-order(correctedFCs[,toPlot[j]])
            predictions1<-rownames(correctedFCs)[O1]   
            XLAB='sgRNA logFC percentile'
        }else{
            sgProf<-correctedFCs[,toPlot[j]]
            names(sgProf)<-rownames(correctedFCs)
            gProf<-ccr.geneMeanFCs(sgProf,libraryAnnotation)
            O1<-order(gProf)
            predictions1<-names(gProf)[O1]
            XLAB='gene Avg logFC percentile'
        }
        
        toTest<-c("BAGEL nonEssential","BAGEL Essential",
                  paste("Amp (PNs >= ",minCN,")",sep=''),
                  paste("Amp (PNs >= ",minCN,") notExp",sep=''))
        
        AUCs<-vector()
        for (i in 1:length(toTest)){
            currentSet<-guideSets[[toTest[i]]]
            currentSet<-intersect(currentSet,predictions1)
            pp1<-cumsum(is.element(predictions1,currentSet))/length(currentSet)
            
            if(i == 1){
                plot(100*(1:length(predictions1))/length(predictions1),
                     100*pp1,type='l',xlim=c(0,100),ylim=c(0,100),ylab='% Recall',
                     col=COLS[i],xlab=XLAB,
                     main=MAIN,lwd=2)
            }else{
                lines(100*(1:length(predictions1))/length(predictions1),
                      100*pp1,type='l',xlim=c(0,100),ylim=c(0,100),col=COLS[i],lwd=2)
            }
            
            AUCs[i]<-trapz(1:length(predictions1)/length(predictions1),pp1)
        }
        
        AUCcombo<-rbind(AUCcombo,AUCs)
        if (j == 1){
            legend('bottomright',toTest,cex=0.8,lwd=2,col=COLS)
        }
    }
    
    rownames(AUCcombo)<-c('pre-CRISPRcleanR','post-CRISPRcleanR')
    colnames(AUCcombo)<-toTest
    
    return(AUCcombo)
}


## other exported non documented functions

#### Assessment and visualisation
ccr.impactOnPhenotype<-function(MO_uncorrectedFile,
                                MO_correctedFile,
                                sigFDR=0.05,
                                expName='expName',
                                display=TRUE){
    
    th<-sigFDR
    
    pre<-read.table(MO_uncorrectedFile,
                    sep='\t',
                    header=TRUE, row.names = 1,stringsAsFactors = FALSE)
    
    pre<-pre[order(rownames(pre)),]
    
    post<-read.table(MO_correctedFile,
                     sep='\t',
                     header=TRUE, row.names = 1,stringsAsFactors = FALSE)
    
    post<-post[order(rownames(post)),]
    
    preD<-which(pre$neg.fdr<sigFDR & pre$pos.fdr>=sigFDR)
    preE<-which(pre$neg.fdr>=sigFDR & pre$pos.fdr<sigFDR)
    preNULL<-setdiff(1:nrow(pre),c(preD,preE))    
    
    postD<-which(post$neg.fdr<sigFDR & post$pos.fdr>=sigFDR)
    postE<-which(post$neg.fdr>=sigFDR & post$pos.fdr<sigFDR)
    postNULL<-setdiff(1:nrow(post),c(postD,postE))    
    
    aDD<-length(intersect(preD,postD))
    aDN<-length(intersect(preD,postNULL))
    aDE<-length(intersect(preD,postE))
    
    aND<-length(intersect(preNULL,postD))
    aNN<-length(intersect(preNULL,postNULL))
    aNE<-length(intersect(preNULL,postE))
    
    aED<-length(intersect(preE,postD))
    aEN<-length(intersect(preE,postNULL))
    aEE<-length(intersect(preE,postE))
    
    cm<-matrix(c(aDD,aDN,aDE,aND,aNN,aNE,aED,aEN,aEE),3,3,dimnames = list(c('cD','cN','cE'),c('uD','uN','uE')))
    cm[is.na(cm)]<-0
    
    IMPACTEDg<-100*sum(triu(cm,1)+tril(cm,-1))/sum(c(cm))
    IMPACTED_phenGenes<-100*(cm[2,1]+cm[2,3]+cm[3,1]+cm[3,2])/sum(c(cm[,c(1,3)]))
    DISTORTEDg<-100*(cm[1,3]+cm[3,1])/sum(c(cm))
    DISTORTED_phenGenes<-100*(cm[1,3]+cm[3,1])/sum(c(cm[,c(1,3)]))
    
    geneCounts<-cm
    
    colnames(cm)<-paste(colSums(cm),c('loss of fitness','no phenotype','gain of fitness'),sep='\n')
    cm<-cm/t(matrix(rep(colSums(cm),nrow(cm)),3,3))
    
    
    if(display){
        par(mar=c(5,4,4,10))
        par(xpd=TRUE)
        barplot(100*cm,col=c('red','gray','blue'),border = FALSE,main=expName,ylab='%',xlab='original counts')
        legend('right',c('loss of fitness','no phenotype','gain of fitness'),inset = c(-.5,0),title = 'Corrected counts',
               fill=c('red','gray','blue'),border = NA)
        
        
        par(mfrow=c(2,2))
        par(mar=c(0,0,2,0))
        pie(c(IMPACTEDg,100-IMPACTEDg),col=c('blue','white'),
            border = 'gray',
            labels = c(paste(format(IMPACTEDg,digits=4),'%',sep=''),''),
            main = 'Overall impact')
        
        pie(c(DISTORTEDg,100-DISTORTEDg),col=c('blue','white'),
            border = 'gray',
            labels = c(paste(format(DISTORTEDg,digits=4),'%',sep=''),''),
            main = 'Overall distortion')
        
        pie(c(IMPACTED_phenGenes,100-IMPACTED_phenGenes),col=c('darkgreen','white'),
            border = 'gray',
            labels = c(paste(format(IMPACTED_phenGenes,digits=4),'%',sep=''),''),
            main = 'Impact (G/L fitness genes)')
        
        pie(c(DISTORTED_phenGenes,100-DISTORTED_phenGenes),col=c('darkgreen','white'),
            border = 'gray',
            labels = c(paste(format(DISTORTED_phenGenes,digits=4),'%',sep=''),''),
            main = 'Distortion (G/L fitness genes)')
        
        
    }
    
    
    dimnames(geneCounts)<-list(`corrected counts`=c('dep.','null','enr.'),`original counts`=c('dep.','null','enr.'))
    
    
    id<-intersect(preD,postE)
    to_bind<-cbind(pre[id,c('neg.fdr','pos.fdr')],post[id,c('neg.fdr','pos.fdr')])
    
    id<-intersect(preE,postD)
    to_bind<-rbind(to_bind,cbind(pre[id,c('neg.fdr','pos.fdr')],post[id,c('neg.fdr','pos.fdr')]))
    colnames(to_bind)<-paste(c('','','ccr.','ccr.'),colnames(to_bind),sep='')
    
    id<-intersect(preD,postD)
    to_bind_c<-cbind(pre[id,c('neg.fdr','pos.fdr')],post[id,c('neg.fdr','pos.fdr')])
    
    id<-intersect(preE,postE)
    to_bind_c<-rbind(to_bind_c,cbind(pre[id,c('neg.fdr','pos.fdr')],post[id,c('neg.fdr','pos.fdr')]))
    colnames(to_bind_c)<-paste(c('','','ccr.','ccr.'),colnames(to_bind_c),sep='')
    
    
    return(list(`GW_impact %`=IMPACTEDg,
                `Phenotype_G_impact %`=IMPACTED_phenGenes,
                `GW_distortion %`=DISTORTEDg,
                `Phenotype_G_distortion %`=DISTORTED_phenGenes,
                geneCounts=geneCounts,
                distortion=to_bind))
}
### Utils


ccr.ExecuteMageck<-function(mgckInputFile,expName='expName',normMethod='none',
                              outputPath='./'){
    fc<-read.table(mgckInputFile,sep='\t',header=TRUE)
    

    Cnames<-colnames(fc)[3]
    Tnames<-colnames(fc)[4:ncol(fc)]
    
    textbunch<-'mageck test -k'
    textbunch<-paste(textbunch,' ',mgckInputFile,' -c ',paste(Cnames,collapse=','),' -t ',paste(Tnames,collapse=','),' -n ',
                     outputPath,expName,' --norm-method ',normMethod,sep='')
    system(textbunch)
    
    geneSummaryFN<-paste(expName,'.gene_summary.txt',sep='')
    
    return(geneSummaryFN)
}


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
ccr.cohens_d <- function(x, y) {
    lx <- length(x)- 1
    ly <- length(y)- 1
    
    md  <- abs(mean(x,na.rm = TRUE) - mean(y,na.rm = TRUE))        ## mean difference (numerator)
    csd <- lx * var(x,na.rm = TRUE) + ly * var(y,na.rm = TRUE)
    
    csd <- csd/(lx + ly)
    csd <- sqrt(csd)                     ## common sd computation
    
    cd  <- md/csd                        ## cohen's d
    return(cd)
}
ccr.get.guideSets<-function(cellLine,GDSC.geneLevCNA=NULL,CCLE.gisticCNA=NULL,RNAseq.fpkms=NULL,libraryAnnotation){
    if(is.null(GDSC.geneLevCNA)){
        data(GDSC.geneLevCNA,envir = environment()) 
    }
    if(is.null(CCLE.gisticCNA)){
        data(CCLE.gisticCNA,envir = environment()) 
    }
    if(is.null(RNAseq.fpkms)){
        data(RNAseq.fpkms,envir = environment()) 
    }
    
    GDSC.cna<-ccr.get.gdsc1000.AMPgenes(cellLine = cellLine,GDSC.geneLevCNA = GDSC.geneLevCNA,minCN = 0)
    CCLE.cna<-ccr.get.CCLEgisticSets(cellLine = cellLine,CCLE.gisticCNA = CCLE.gisticCNA)
    notExp<-ccr.get.nonExpGenes(cellLine = cellLine,RNAseq.fpkms = RNAseq.fpkms)
    
    data(BAGEL_essential,envir = environment()) 
    data(BAGEL_nonEssential,envir = environment()) 
    
    data(EssGenes.DNA_REPLICATION_cons,envir = environment()) 
    data(EssGenes.KEGG_rna_polymerase,envir = environment()) 
    data(EssGenes.PROTEASOME_cons,envir = environment()) 
    data(EssGenes.ribosomalProteins,envir = environment()) 
    data(EssGenes.SPLICEOSOME_cons,envir = environment()) 
    
    CFEgenes<-unique(c(EssGenes.DNA_REPLICATION_cons,EssGenes.KEGG_rna_polymerase,
                       EssGenes.PROTEASOME_cons,EssGenes.ribosomalProteins,EssGenes.SPLICEOSOME_cons))
    
    BAGEL_essOnly<-setdiff(BAGEL_essential,CFEgenes)
    
    options(warn=-1)
    BAGEL_essential<-ccr.genes2sgRNAs(BAGEL_essential,libraryAnnotation = libraryAnnotation)
    BAGEL_nonEssential<-ccr.genes2sgRNAs(BAGEL_nonEssential,libraryAnnotation = libraryAnnotation)
    EssGenes.DNA_REPLICATION_cons<-ccr.genes2sgRNAs(EssGenes.DNA_REPLICATION_cons,libraryAnnotation = libraryAnnotation)
    EssGenes.KEGG_rna_polymerase<-ccr.genes2sgRNAs(EssGenes.KEGG_rna_polymerase,libraryAnnotation = libraryAnnotation)
    EssGenes.PROTEASOME_cons<-ccr.genes2sgRNAs(EssGenes.PROTEASOME_cons,libraryAnnotation = libraryAnnotation)
    EssGenes.ribosomalProteins<-ccr.genes2sgRNAs(EssGenes.ribosomalProteins,libraryAnnotation = libraryAnnotation)
    EssGenes.SPLICEOSOME_cons<-ccr.genes2sgRNAs(EssGenes.SPLICEOSOME_cons,libraryAnnotation = libraryAnnotation)
    
    CFEgenes<-ccr.genes2sgRNAs(CFEgenes,libraryAnnotation = libraryAnnotation)
    BAGEL_essOnly<-ccr.genes2sgRNAs(BAGEL_essOnly,libraryAnnotation = libraryAnnotation)
    
    DelGDSC<-ccr.genes2sgRNAs(libraryAnnotation,GDSC.cna$Gene[GDSC.cna$CN==0])
    DelCCLE<-ccr.genes2sgRNAs(libraryAnnotation,CCLE.cna$gm2)
    notExpG<-ccr.genes2sgRNAs(libraryAnnotation,notExp)
    gistic1<-ccr.genes2sgRNAs(libraryAnnotation,CCLE.cna$gp1)
    gistic2<-ccr.genes2sgRNAs(libraryAnnotation,CCLE.cna$gp2)
    gistic1ne<-ccr.genes2sgRNAs(libraryAnnotation,intersect(CCLE.cna$gp1,notExp))
    gistic2ne<-ccr.genes2sgRNAs(libraryAnnotation,intersect(CCLE.cna$gp2,notExp))
    
    pn2<-ccr.genes2sgRNAs(libraryAnnotation,GDSC.cna$Gene[GDSC.cna$CN>=2])
    pn4<-ccr.genes2sgRNAs(libraryAnnotation,GDSC.cna$Gene[GDSC.cna$CN>=4])
    pn8<-ccr.genes2sgRNAs(libraryAnnotation,GDSC.cna$Gene[GDSC.cna$CN>=8])
    pn10<-ccr.genes2sgRNAs(libraryAnnotation,GDSC.cna$Gene[GDSC.cna$CN>=10])
    
    pn2ne<-ccr.genes2sgRNAs(libraryAnnotation,intersect(GDSC.cna$Gene[GDSC.cna$CN>=2],notExp))
    pn4ne<-ccr.genes2sgRNAs(libraryAnnotation,intersect(GDSC.cna$Gene[GDSC.cna$CN>=4],notExp))
    pn8ne<-ccr.genes2sgRNAs(libraryAnnotation,intersect(GDSC.cna$Gene[GDSC.cna$CN>=8],notExp))
    pn10ne<-ccr.genes2sgRNAs(libraryAnnotation,intersect(GDSC.cna$Gene[GDSC.cna$CN>=10],notExp))
    
    options(warn=0)
    
    guideSets<-list(DelGDSC,
                    DelCCLE,
                    notExpG,
                    gistic1,
                    gistic2,
                    gistic1ne,
                    gistic2ne,
                    pn2,
                    pn4,
                    pn8,
                    pn10,
                    pn2ne,
                    pn4ne,
                    pn8ne,
                    pn10ne,
                    CFEgenes,
                    BAGEL_essential,
                    BAGEL_essOnly,
                    BAGEL_nonEssential)
    
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
              'BAGEL nonEssential')
    names(guideSets)<-nnames
    return(guideSets)
}
ccr.fixedFDRthreshold<-function(FCsprofile,TruePositives,TrueNegatives,th){
    presentGenes<-intersect(c(TruePositives,TrueNegatives),names(FCsprofile))
    predictions<-FCsprofile[presentGenes]
    observations<-is.element(presentGenes,TruePositives)+0
    names(observations)<-presentGenes
    RES<-roc(observations,predictions,'>')
    COORS<-coords(RES,'all',ret = c('threshold','ppv'))
    FDRpercTh<-max(COORS['threshold',which(COORS['ppv',]>=(1-th))])
    FDRpercRANK<-max(which(sort(FCsprofile)<=FDRpercTh))
    
    return(list(FCth=FDRpercTh,RANK=FDRpercRANK))
}
ccr.sd_guideFCs<-function(FCs,distorted_genes,libraryAnnotation){
    
    SDS<-aggregate(FCs,by=list(libraryAnnotation[names(FCs),'GENES']),'sd')
    XSDS<-SDS$x
    names(XSDS)<-SDS$Group.1
    
    boxplot(XSDS[setdiff(names(XSDS),distorted_genes)],XSDS[distorted_genes])
    
    print(t.test(XSDS[setdiff(names(XSDS),distorted_genes)],XSDS[distorted_genes]))
    
    
}


