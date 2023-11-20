library(data.table)
library(vcfR)
library(dplyr)

id <- read.csv("nwdid.csv")
dat_vcf <- fread("mutect2.mutect2.05.merge.vcf.gz")

## extract samples
col_id <- c(colnames(dat_vcf)[1:9], colnames(dat_vcf)[colnames(dat_vcf) %in% id$NWD_ID])
dat <- dat_vcf %>% select(all_of(col_id))

# reference sequence
Ref <- read.csv("mt_reference.csv",stringsAsFactors = F)
.mtRef <- Ref$ref
.mtLength <- nrow(Ref)

# remove insertion and deletion
dat <- dat[nchar(dat$ALT) == nchar(dat$REF),] #23323 --> 20912

# remove POS with reference allele = alternative allele
dat <- dat[dat$REF != dat$ALT,] #20912 --> 16334

# only maintain same reference allele as the reference sequence
dat2 <- data.frame(matrix(ncol = ncol(dat), nrow = 0))
colnames(dat2) <- colnames(dat)
for (i in unique(dat$POS)){
  a <- dat[dat$POS==i,]
  a <- a[a$REF %in% Ref$ref[i],]
  dat2 <- rbind(dat2,a)
} #16334 --> 16236


#### create allele, coverage, frequency matrixes #### 
mtReadVCF_raw <- function(VCF_file, n_sample){
  dat_vcf <- VCF_file[,2:ncol(VCF_file)]
  pos_var <- unique(dat_vcf$POS)
  pos_nonvar <- setdiff(c(1:.mtLength), dat_vcf$POS) #find all rows in ref that aren't in dat_vcf
  allele_pos_nonvar <- .mtRef[pos_nonvar]
  allele_pos_nonvar_mat <-replicate(n_sample, allele_pos_nonvar)

  allele <- matrix("0", nrow=.mtLength, ncol=n_sample)
  freq <- matrix("1", nrow=.mtLength, ncol=n_sample)
  coverage <- matrix(NA, nrow=.mtLength, ncol=n_sample)
  rownames(allele) <- seq_len(.mtLength)
  rownames(freq) <- seq_len(.mtLength)
  rownames(coverage) <- seq_len(.mtLength)
  subjects_index <- c((ncol(dat_vcf)-n_sample+1):ncol(dat_vcf))
  colnames(allele) <- colnames(dat_vcf)[subjects_index]
  colnames(freq) <- colnames(dat_vcf)[subjects_index]
  colnames(coverage) <- colnames(dat_vcf)[subjects_index]

  allele[pos_nonvar,] <- allele_pos_nonvar_mat #fill rows only in ref with ref allele

  # split the values of allele, freq and coverage from VCF
  dat_vcf_geno <- dat_vcf[,(ncol(dat_vcf)-n_sample+1):ncol(dat_vcf)]

  for(j in pos_var){
    # for each position in VCF file
    a <- dat_vcf[dat_vcf$POS==j,]
    allele_a <- matrix(nrow = nrow(a), ncol = n_sample)
    freq_a <- matrix(nrow=nrow(a), ncol=n_sample)
    coverage_a <- matrix(nrow=nrow(a), ncol=n_sample)
    colnames(allele_a) <- colnames(dat_vcf)[subjects_index]
    colnames(freq_a) <- colnames(dat_vcf)[subjects_index]
    colnames(coverage_a) <- colnames(dat_vcf)[subjects_index]
    b <- a[,(ncol(a)-n_sample+1):ncol(a)]

    for (k in 1:nrow(a)){
      pos_alleles <- c(a$REF[k], unlist(strsplit(a$ALT[k], split="")))
      x <- sapply(strsplit(as.character(b[k,]), ":"), "[", c(1,2,3))

      # allele
      alleles_point <- x[1, ]
      alleles_point[alleles_point=="."] <- ""
      ##cannot use "|" here because R understand "|" as or, then R will split string with ""
      alleles_split <- ifelse(grepl("/", alleles_point), strsplit(alleles_point, split = "/"), strsplit(alleles_point, split = "\\|"))
      alleles_split <- lapply(alleles_split, unique)
      alleles_split <- lapply(alleles_split, as.numeric)
      alleles_split <- sapply(alleles_split, function(x) x<-x+1)
      #make sure reference allele is the first allele, because 0|1 or 1|0 are different (phased genotype) (0/1, 1/0 are unphased genotype)
      alleles_split <- sapply(alleles_split, sort)

      for (i in seq_len(n_sample)){
        alleles_split_letter <- paste(pos_alleles[alleles_split[[i]]], collapse = "/")
        allele_a[k, i] <- alleles_split_letter
        #freq (if genotype is NA, freq and coverage are also NA)
        freq_a[k, i] <- ifelse(x[3,i]==".", NA, ifelse(x[3,i]==1, 1, as.numeric(x[3, i])))
        #coverage
        coverage_a[k,i] <- ifelse(x[2,i]==".", NA, as.numeric(x[2,i]))
      }
    }

    ## put allele_a, seq_a, coverage_p to allele, seq, coverage (arrange one or multiple rows for each gene to one row)
    allele_a[allele_a == ""] <- NA
    for (p in seq_len(n_sample)){
      if (colSums(is.na(allele_a))[p] == nrow(a)){
        allele[j,p] <- NA
        freq[j,p] <- NA
        coverage[j,p] <- NA
      }else{
        # some sample may have >1 alternative alleles
        allele[j,p] <- paste(unique(unlist(strsplit(paste0(na.omit(allele_a[,p]), collapse = "/"), split="/"))), collapse = "/")
        freq[j,p] <- paste(c(1-sum(as.numeric(freq_a[,p]), na.rm = T), paste0(na.omit(freq_a[,p]))), collapse = "/")
        coverage[j,p] <- unique(coverage_a[!is.na(coverage_a[,p])][p])
      }
    }
  }

  # fill NA in freq with 1, and allele with reference allele in the ref dataset
  freq[is.na(freq)] <- 1
  for (z in pos_var){
    allele[z,is.na(allele[z,])] <- Ref[z,"ref"]
  }
  # output
  out_list <- list("allele" = allele, "freq" = freq, "coverage" = coverage)
  return(out_list)
}


vcf_raw <- mtReadVCF_raw(VCF_file=dat2, n_sample=n) #n = number of samples in dat2


save(vcf_raw, file = "extract_vcf.RData")



