######################################
#### calcuate heteroplasmy burden ####
######################################
library(dplyr)



## read in data
mito_dat <- read.table("MLC_score.tsv")
load("aaf_fhs.RData")


rownames(AAF.m) <- rownames(allele)
aaf <- AAF.m
#allele[1:10, 1:10]
#freq[1:10, 1:10]
#coverage[1:10, 1:10]
#aaf[1:10, 1:10]



## QC: set the AAF with coverage<250 to 0
#check <- coverage[coverage<250 & is.na(coverage)==F]
aaf[coverage<250 & is.na(coverage)==F] <- NA



# not reliable loci (NUMT) --> remove
loci_NUMT <- c(1:61,301,302,310,316,499,567,3107,16088:16569)
aaf <- aaf[-loci_NUMT, ] # 16569 --> 16019, 7210 samples
coverage <- coverage[-loci_NUMT, ]
freq <- freq[-loci_NUMT, ]
allele <- allele[-loci_NUMT, ]


## calculated mito-score weighted heteroplasmy
# actually there is no need to set threshold 0.05-0.95, since mitoHPC has already applied it, 
# but we still need to identify those variants with AAF with 0 1
thre_spec_lower <- 0.05
thre_spec_upper <- 0.95
allele <- allele[,colnames(aaf)]
freq <- freq[,colnames(aaf)]
aaf_cat0 <- ((aaf >= thre_spec_lower) & (aaf <= thre_spec_upper))
#aaf_cat0[1:10, 1:10]
aaf_cat0 <- aaf_cat0 + 0 #transfer FALSE to 0
#aaf_cat0[1:10, 1:10]
#dim(aaf_cat0)

aaf_cat0_nomiss <- aaf_cat0
#check2 <- aaf_cat0_nomiss[is.na(aaf_cat0_nomiss)]
aaf_cat0_nomiss[is.na(aaf_cat0_nomiss)] <- 0
sit_aaf_cat0 <- rownames(aaf_cat0)

heter_count_loci <- rowSums(aaf_cat0, na.rm = TRUE)
sit_loop_check <- as.numeric(sit_aaf_cat0[heter_count_loci>0]) #pos with heteroplasmy

aaf_cat1 <- aaf_cat0
pos_mito_miss <- setdiff(x=c(1:16569), y=unique(mito_dat$POS))
#check3 <- unique(mito_dat$POS)
sit_loop_check <- sit_loop_check[!sit_loop_check%in%pos_mito_miss] #position with heteroplasmy missing in mito_dat

#check4 <- seq_along(sit_loop_check) #1:1652
for (i in seq_along(sit_loop_check)) { #row (pos)
  site <- sit_loop_check[i]
  row_index_aaf <- match(as.character(site), rownames(aaf_cat0)) # the rownames(aaf_cat0)[row_index_aaf] == sit_loop_check[i]
  
  #check5 <- colSums(aaf_cat0, na.rm = TRUE)
  for(j in 1:ncol(aaf)) { #col (num indivi)
    
    if(aaf_cat0_nomiss[row_index_aaf,j] ==1) { #if all variant #HOW ABOUT IF NOT EQUAL TO 1??
      
      alt <- allele[row_index_aaf,j]
      
      if(nchar(alt)==1){
        alt2 <- alt
      }else if (nchar(alt)>1) { #if >1 alleles present, find the allele corresponding to aaf
        all_alts  <- strsplit(alt,"/")[[1]]
        
        all_freqs <- as.numeric(strsplit(freq[row_index_aaf,j],"/")[[1]])
        alt_freq  <- aaf[row_index_aaf,j]
        alt_index <- which(all_freqs==alt_freq)
        
        ref <- mito_dat[mito_dat$POS==site,]$REF[1]
        alt2 <- all_alts[alt_index]
        
        if (length(alt2)>=2){
          alt2 <- alt2[!alt2 %in% c(ref,"I","D")]
          alt2 <- alt2[1]
        }
        
        if(length(alt2)!=1 | length(all_freqs)!=length(all_alts)) {
          stop(paste0("r:",site,";c:",j,";length(alt2)!=1",alt2))
        }
      }
      
      mitoscore <- mito_dat[mito_dat$POS==site & mito_dat$ALT==alt2,]$mito_lc_score
      aaf_cat1[row_index_aaf,j] <- mitoscore
    }
  }
}



## heteroplasmic burden
weight_heter <- colSums(aaf_cat1, na.rm = T)  #sum entries in each column separately, remove missing val
weight_heter <- as.data.frame(weight_heter)
weight_heter$NWD <- rownames(weight_heter)

raw_heter <- colSums(aaf_cat0_nomiss, na.rm = T)
raw_heter <- as.data.frame(raw_heter)
raw_heter$NWD <- rownames(raw_heter)

hetero_burden <- left_join(weight_heter, raw_heter, by = "NWD")
hetero_burden$bi_heter <- ifelse(hetero_burden$raw_heter==0, "0", "1")
hetero_burden <- hetero_burden[,c("NWD", "raw_heter", "bi_heter", "weight_heter")]


#### save ####
write.csv(hetero_burden, file = "het_score.csv", row.names = F)


#### SAMPLE REMOVAL ####
# I haven't included sample removal here. after getting het_score.csv, we can remove samples with mtDNA CN < 40; contamination level > 3%; raw_heter > 5;

