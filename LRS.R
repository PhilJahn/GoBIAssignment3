LRS <- function(incl, total, group){
  
  coln <- length(group)
  
  inc0 <- sum(incl)
  tot0 <- sum(total)
  
  inc1 <- 0
  tot1 <- 0
  
  inc2 <- 0
  tot2 <- 0
  
  for(i in 1:coln){
    inci <- incl[i]
    toti <- total[i]
    
    if(group[i] == 1){
      inc1 <- inc1 + inci
      tot1 <- tot1 + toti
    }
    else{
      inc2 <- inc2 + inci
      tot2 <- tot2 + toti    
    }

  }
  
  p0 <- inc0/tot0
  p1 <- inc1/tot1
  p2 <- inc2/tot2
  
  llreduced <- 0
  llfull <- 0
  
  for(j in 1:coln){
    incj <- incl[j]
    totj <- total[j]
    llreduced <- llreduced + log(dbinom(incj,totj,p0))
    
    if(group[j] == 1){
      llfull <- llfull + log(dbinom(incj,totj,p1))
    }
    else{
      llfull <- llfull + log(dbinom(incj,totj,p2))   
    }
  }
  
  lrs <- -2 * (llreduced - llfull)
  
  pvalue <- pchisq(lrs, df=1, lower.tail=FALSE)
  
  result <- list(p0 = p0, p1 = p1, p2=p2, llreduced = llreduced, llfull = llfull, lrs = lrs, pvalue = pvalue)
  return(result)
}

diff.splicing <- function(psi.files, group){
  
  filecount <- length(psi.files)

  file <- read.table(psi.files[1],header=TRUE, sep="\t",stringsAsFactors=FALSE)
  
  content <- data.frame(gene = file$gene,exon = file$exon, incl1= file$num_incl_reads, total1 = file$num_total_reads)
  
  for(i in 2:filecount){
    file <- read.table(psi.files[i],header=TRUE, sep="\t")

    
    newcontent <- data.frame(gene = file$gene, exon =file$exon, file$num_incl_reads, file$num_total_reads)
    colnames(newcontent)[3] <- paste("incl",i, sep="")
    colnames(newcontent)[4] <- paste("total",i, sep="")
    
    content <- merge(content, newcontent, by=c("gene","exon"), all=TRUE)
  }
  
  content[is.na(content)] <- 0
  
  exonCount <- nrow(content)
  
  incl <- content[seq(3, ncol(content), 2)]
  
  incl <- sapply(incl, as.numeric)
  
  total <- content[seq(4, ncol(content), 2)]
  total <- sapply(total, as.numeric)
  
  lrs <- data.frame(LRS(incl[1,],total[1,],group))
  
  for(i in 2:exonCount){
    lrs <- rbind(lrs,LRS(incl[i,],total[i,],group))
  }
  
  output <- data.frame(content[1],content[2])
  
  output <- cbind(output,lrs)
  
  output <- cbind(output, padj=p.adjust(output$pvalue, method = "BH"))

  return(output)
}