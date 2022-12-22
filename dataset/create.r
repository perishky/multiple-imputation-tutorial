source("geo.r")

#################################################

## dataset 1 samples
samples <- geo.samples("GSE50660")
samples <- with(samples, {
    data.frame(gsm=geo_accession,
               gse="GSE50660",
               smoking=get.characteristic(characteristics_ch1,"smoking"),
               age=get.characteristic(characteristics_ch1,"age","numeric"),
               sex=get.characteristic(characteristics_ch1,"gender"),
               stringsAsFactors=F)
})
samples$smoking <- recode(samples$smoking,
                          c("0"="never",
                            "1"="former",
                            "2"="current"))
samples$sex <- tolower(samples$sex)
rownames(samples) <- NULL

## dataset 1 methylation
link <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50660/suppl/GSE50660_matrix_processed.txt.gz"
filename <- basename(link)
if (!file.exists(filename))
    download.file(link, filename)
meth <- read.table(filename, header=T,sep="\t",quote="",comment.char="") ## 5 minutes

rownames(meth) <- meth$ID_REF
meth <- as.matrix(meth[,-1])

## dataset 1 link samples and methylation
colnames(meth) <- samples$gsm

######################################

library(meffil) ## https://github.com/perishky/meffil

reference <- "blood gse35069 complete"
cc <- meffil.estimate.cell.counts.from.betas(meth, reference, verbose=T)
## 2 minutes

## sva for age
mod <- model.matrix(~ age + smoking + sex, samples)
mod <- cbind(mod, cc)
set.seed(20171025)
system.time(sva.age <- sva(meth, mod = mod, mod0=subset(mod, select=-age)))
## 5 hours
sva.age <- sva.age$sv
rownames(sva.age) <- colnames(meth)

## sva for current smoking
idx <- which(samples$smoking != "former")
mod <- model.matrix(~ age + smoking + sex, samples[idx,])
mod <- cbind(mod, cc[idx,])
mod[,"smokingnever"] <- 1-mod[,"smokingnever"]
colnames(mod)[match("smokingnever", colnames(mod))] <- "smoking"
set.seed(20171025)
system.time(sva.smoking <- sva(meth[,idx], mod=mod, subset(mod, select=-smoking)))
## 1 hour
sva.smoking <- sva.smoking$sv
rownames(sva.smoking) <- colnames(meth)[idx]

sites <- readLines("sites.txt")
meth <- meth[sites,]

samples$smoking <- as.factor(samples$smoking)

## save dataset 1
save(samples, meth, cc, sva.age, sva.smoking, file="dataset.rda")


