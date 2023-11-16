load("extract_vcf.RData")

allele <- vcf_raw$allele
freq <- vcf_raw$freq
coverage <- vcf_raw$coverage

Ref <- read.csv("mt_reference.csv", stringsAsFactors = F)

### calculate AAF ####
ref2 <- Ref[,2]
ref2 <- as.character(ref2)
subjectID <- colnames(allele)

# make sure allele, coverage, and freq have same order of people
freq <- freq[,colnames(allele)]
coverage <- coverage[,colnames(allele)]

AAF.m <- array(0, dim(allele))
system.time({
  for (m in 1:dim(allele)[2]) {
    
    m.allele   <- allele[,m]
    m.freq     <- freq[,m]
    
    complex.allele <- grep("/", m.allele) #return row numbers for rows containing "/"
    AAF <- array(NA, 16569) #16569: number of genes
    # Xianbang: I added a if statement here since complex.allele may have length of 0, it is very rare in practice, but is possible theoretically
    if(length(complex.allele)==0){
      AAF <- ifelse(m.allele == ref2, 0, 1)
    }else{
      AAF[-complex.allele] <- ifelse(m.allele[-complex.allele] == ref2[-complex.allele], 0, 1)
      
      complex.all <- strsplit(m.allele[complex.allele],split="/")
      complex.freqsplit <- strsplit(m.freq[complex.allele],split="/")
      complex.freqsplit <- lapply(complex.freqsplit, as.numeric)
      complex.ref <- ref2[complex.allele]
      AAF[complex.allele] <- mapply(function(x, y, z){ pos <- which(x %in% z); if(length(pos) > 0) max(y[-pos]) else max(y)  } , x= complex.all, y = complex.freqsplit, z = complex.ref)
    }    #if samples at specific position contain reference allele, length(pos)=1, if no reference allele (different alternative allele), length(pos)=0
    
    AAF.m[,m]<-AAF
  }
  
  # Katia: do we really need this?
  colnames(AAF.m) <- subjectID
  # Xianbang: we can assign the ID to the colnames of PAA matrix
  
})

save(allele, freq, coverage, AAF.m, file = "aaf_fhs.RData")
