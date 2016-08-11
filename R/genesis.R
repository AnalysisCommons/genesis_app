#== Args 
args<-commandArgs(TRUE)
#===mandatory parameters
phenotype.file <- args[1]
outcome.name <- args[2]
outcome.type <-  args[3]
covariate.string <- args[4]
snpinfo.file <- args[5]
genotype.files <- args[6]
output.file <- args[7]

#==optional parameters
kinship.matrix <- args[8]
gender.column <- args[9]
pheno.id <- args[10]
nsmatch <- args[11]


# added these to JSON
BUFFER <- as.numeric(args[12]) #jb
gene.file <- args[13] #jb
snp.filter <- args[14] #jb
gene.filter <- args[15] #jb
top.maf <- as.numeric(args[16]) #jb
test.requested <-  args[17]
burden.test.requested <-  args[18]
min.mac <- as.integer(args[19])

# GLOBAL VARIABLES
collapsing.tests <- c("SKAT", "Burden")
burden.tests <- c("Score", "Wald", "Firth")
single.variant.tests <- c("Wald", "Score")
GetFamilyDistribution <- function(response.type) {
               if (response.type == "Continuous"){
                      family = "gaussian"
               } else if (response.type == "Dichotomous"){
                      family = "binomial"
               } else {
                      msg = paste("Don't know how to deal with response type", response.type)
                      stop(msg)
               }
               return(family)
           }

getMAC <- function(gds.file){
  geno.dat <- seqGetData(gds.file,"genotype")
  apply(geno.dat,3,function(x){min(sum(x == 0,na.rm=T),sum(x == 1,na.rm=T))})
}

cat('output.file',output.file,'\n')
cat('kinship.matrix',kinship.matrix,'\n')
cat('gender.column',gender.column,'\n')
cat('nsmatch',nsmatch,'\n')
cat('buffer',BUFFER,'\n')
cat('gene.file',gene.file,'\n')
cat('snp.filter',snp.filter,'\n')
cat('top.maf',top.maf,'\n')
cat('test.requested',test.requested,'\n')
cat('outcome.type',outcome.type,'\n')

if (!(test.requested %in% c(collapsing.tests,single.variant.tests))){
     msg = paste("The requested test:", test.requested, "is not available!")
     stop(msg)
}


library(SeqArray)
library(SeqVarTools)
library(GWASTools)
library(gap)
library(Matrix)
library(plyr)
library(gdsfmt)
library(bdsmatrix)
library(parallel)
library(GENESIS)
#setMKLthreads(1)
library(data.table)
library(doMC)

num.cores <- detectCores(logical=TRUE)
registerDoMC(cores=num.cores-1)
cat('Number of cores', num.cores,'\n')

## Setup
source("/tmp/pipelineFunctions.R")
covariates <- split.by.comma(covariate.string)  

##  Gene list
cat('reading gene file...')

cat('Timings with read.table....')
system.time({ kg = read.table(gene.file, as.is=T, sep=',', header=T) })
cat('Timings with fread....')
system.time({ kg = fread(gene.file, stringsAsFactors=F, sep=',', header=T) })
cat('GENE Filter',gene.filter,'\n')
kg = eval(parse(text= paste0('subset(kg,',gene.filter,')')))

cat(NROW(kg),'done\n')
genes <- kg$name


## snp info
cat('Reading snpinfo....')
snpinfo <- fread(snpinfo.file)
cat('done\n')

## Filtering snp info
## Ive heard this is a horrible was to program things... but I do this sometimes -- feel free to edit
cat('Input SNPINFO N=',nrow(snpinfo),'\n')
snpinfo = eval(parse(text= paste0('subset(snpinfo,',snp.filter,')')))
cat('Output SNPINFO N=',nrow(snpinfo),'\n')


## phenotype 
phenotype.data <- read.csv(phenotype.file, header=TRUE, as.is=TRUE)

# Set gender.col to "sex" in both covariates and phenotype
if (gender.column %in% covariates){
    covariates[covariates == gender.column] = "sex"
}
if (gender.column %in% colnames(phenotype.data)){
    colnames(phenotype.data)[colnames(phenotype.data) == gender.column] <- "sex"
}


cat('Input pheno N=',nrow(phenotype.data),'\n')
pheno <- reducePheno(phenotype.data, outcome.name, covariates, pheno.id, "sex")
cat('Output pheno N=',nrow(pheno),'\n')

## Report dropped individuals
dropped.ids.selector <- !(phenotype.data[[pheno.id]] %in% row.names(pheno))
dropped.ids <- phenotype.data[[pheno.id]][dropped.ids.selector] 
if (NROW(dropped.ids) != 0 ) {
    cat("The following ids were dropped because of incomplete cases:", dropped.ids)
}


# For GDS files
f <- seqOpen(genotype.files)
sample.ids <- seqGetData(f, "sample.id")
pheno <- pheno[row.names(pheno) %in% sample.ids,c(outcome.name, covariates, "sex"),drop=F]
full.sample.ids <- sample.ids 

#subset to phenotyped samples
seqSetFilter(f,sample.id = row.names(pheno))

# order pheno to the GDS subject order
sample.ids <- seqGetData(f, "sample.id")
pheno <- pheno[match(sample.ids,row.names(pheno)),,drop=F]



# get position list before any variant filters are applied
pos = seqGetData(f, "position")


## Load KINSHIP matrix
## Kinship doesn't contain all samples
kmatr = get(load(kinship.matrix))
pheno = pheno[row.names(pheno) %in% row.names(kmatr),,drop=F]
kmatr = kmatr[row.names(kmatr) %in% row.names(pheno),colnames(kmatr) %in% row.names(pheno)]
cat('Output pheno in Kinship N=',nrow(pheno),'\n')
kmatr = kmatr[match(row.names(pheno),row.names(kmatr)),match(row.names(pheno),colnames(kmatr))]

# Get sample ids to check order 
seqSetFilter(f,sample.id = row.names(pheno))
sample.ids <- seqGetData(f, "sample.id")

if (!(identical(sample.ids,row.names(pheno)) && identical(row.names(kmatr),row.names(pheno)))){
        stop("Something is off problem with re-ordering")
}

seqClose(f)

sample.data <- data.frame(scanID = row.names(pheno),  
                    pheno, 
                    stringsAsFactors=F)
scan.annotated.frame <- ScanAnnotationDataFrame(sample.data)
modified.pheno = pheno[full.sample.ids,]
row.names(modified.pheno) <- full.sample.ids

sample.data.for.annotated <- data.frame(sample.id = full.sample.ids,
                                        modified.pheno,
                                        stringsAsFactors=F)

annotated.frame <- AnnotatedDataFrame(sample.data.for.annotated)

###################
## NULL MODEL
##################
# Should depend on response type

cat('start fit....\n')
kmatr.ns = as.matrix(kmatr)
nullmod <- fitNullMM(scanData = scan.annotated.frame,
                     outcome = outcome.name,
                     covars = covariates,
                     family = GetFamilyDistribution(outcome.type),
                     covMatList = kmatr.ns)




## For each aggregation unit - named 'gene' in code, but can be any start-stop region
## Work for each gene is parallelized over the cores
mcoptions <- list(preschedule=FALSE, set.seed=FALSE)
sm_obj <- 
foreach (current.gene=genes, 
         .combine=function(...){rbindlist(list(...),fill=T)},
         .multicombine = TRUE,
         .inorder=FALSE,  
         .options.multicore=mcoptions) %dopar% {

  ##############
  ## Apply variant filters to get a list of variants in each gene
  ###############
  
  gidx <- which(kg$name == current.gene)
  if (length(gidx) != 1) {
     msg=paste("The assumptions that the list of genes is unique is violated\n",
               "The gene:", cgene, "repeats", length(gidx), "times")
     stop(msg)
  } 
  
  geneSNPinfo = subset(snpinfo, (POS > (kg[gidx,]$start - BUFFER) & POS < (kg[gidx,]$stop + BUFFER)))
  ## These nsmatch variables are too analysis specific and we may need to move to a array of variants per aggregation unit method
  ## NS=1 Matches nonsyn snps *only* if they match the named gene
  if(nsmatch == '1'){
    geneSNPinfo = subset(geneSNPinfo,  ! ( gene != kg[gidx,]$gene & sc_nonsynSplice == 1))
  ## NS=2 Matches nonsyn snps *only* if they match the named gene AND other variants *only* if they match 'gene2' column in annotation
  }else if(nsmatch == '2'){
    geneSNPinfo = subset(geneSNPinfo,  !( gene != kg[gidx,]$gene & sc_nonsynSplice == 1) || !( gene2 != kg[gidx,]$gene & isCorrDHS == 1) )
  }

  
  snp_idx = which(pos %in%  geneSNPinfo$POS)
  cat(current.gene,'NEW \t',length(snp_idx),'\n')

  #res = list()
  #res[[current.gene]]=list(generes=data.frame())
  if(length(snp_idx) > 0){
    
    ## extract genotypes
    f <- seqOpen(genotype.files)
    
    #seqSetFilter(f,sample.id = row.names(pheno),verbose=FALSE)
    seqSetFilter(f,variant.id = snp_idx,verbose=FALSE)
    
    ## filter to maf and remove monomorphic
    maf <- SeqVarTools::alleleFrequency(f)
    maf <- ifelse(maf < 0.5, maf, 1-maf)
    filtered.alleles = TRUE
    if (test.requested %in% collapsing.tests){
        filtered.alleles <- maf < top.maf
    }
    if (test.requested %in% single.variant.tests){
      mac <- getMAC(f)
      filtered.alleles <- mac  > min.mac
      
    }
    
    seqSetFilter(f,variant.sel=snp_idx[filtered.alleles & maf > 0], verbose=TRUE)
    num.polymorphic.snps <- seqSummary(f, "genotype", check="none", verbose=FALSE)[["seldim"]][3]
    num.snps = length(maf)
    
    cat(num.polymorphic.snps)
    cat('+')
    cat(test.requested)
    
    if(num.polymorphic.snps > 0){
      
       genotype.data <- SeqVarData(f, sampleData=annotated.frame)
       # Collapse test
       if (test.requested %in% collapsing.tests) {
            xlist = list()
            xlist[[1]] = data.frame('variant.id'=seqGetData(f,"variant.id"),
                                    'allele.index'=rep(1,length(seqGetData(f,"variant.id")))
                                    )
            system.time({
                 collapse.results <- assocTestSeq( genotype.data, 
                                         nullmod, 
                                         xlist, 
                                         test=test.requested, 
                                         burden.test=burden.test.requested)
            })
          generes <- cbind(data.frame(gene=current.gene, num.variants= nrow(xlist[[1]])), collapse.results$results)

       } else {  # Single variant test
 
          system.time({generes <- assocTestMM(genotype.data, 
                                              nullmod, 
                                              test = test.requested)})
          generes$gene <- current.gene
          generes$pos=pos[snp_idx[filtered.alleles & maf > 0]] 
      }
   } 
   seqClose(f)
  }
     
  if(!exists("generes")){
    generes= data.frame('gene'=current.gene);
  }
  
  #message(paste(current.gene,"SNPS:",length(snp_idx),"pct:",round(which(genes == current.gene)/length(genes),2),"\t",generes$gene[1],"finished."))
  generes
}
cat("\nFinished Loop\n")
cat("Total output lines: ",nrow(sm_obj),"\n")

  
out <- gzfile(output.file, "w")
write.csv(sm_obj, out, row.names=F)
close(out)
