library(tidyverse)
tx <-as.numeric(Sys.getenv("SGE_TASK_ID"))

gc.correct <- function(coverage, bias) {
    i <- seq(min(bias, na.rm=TRUE), max(bias, na.rm=TRUE), by = 0.001)
    coverage.trend <- loess(coverage ~ bias)
    coverage.model <- loess(predict(coverage.trend, i) ~ i)
    coverage.pred <- predict(coverage.model, bias)
    coverage.corrected <- coverage - coverage.pred + median(coverage)
}

fragpath <- "../fragments"
fragfiles <- list.files(fragpath, pattern=".rds",full.name=TRUE)
fragfile <- fragfiles[tx]

id <- strsplit(basename(fragfile), "\\.")[[1]][1]

outdir <- "." ####
filename <- file.path(outdir, paste0(id, "_bin_100kb.rds"))
if(file.exists(filename)) q('no')

library(GenomicRanges)
library(rtracklayer)
library(Homo.sapiens)
library(BSgenome.Hsapiens.UCSC.hg19)
library(Rsamtools)
class(Homo.sapiens)
library(devtools)
library(biovizBase)
load("./filters.hg19.rda")

library(RCurl)
ABurl <- getURL('https://raw.githubusercontent.com/Jfortin1/HiC_AB_Compartments/master/data/hic_compartments_100kb_ebv_2014.txt', ssl.verifyhost=FALSE, ssl.verifypeer=FALSE)

AB <- read.table(textConnection(ABurl), header=TRUE)
AB <- makeGRangesFromDataFrame(AB, keep.extra.columns=TRUE)

chromosomes <- GRanges(paste0("chr", 1:22),
                       IRanges(0, seqlengths(Hsapiens)[1:22]))

tcmeres <- gaps.hg19[grepl("centromere|telomere", gaps.hg19$type)]##gaps.hg19无法访问这个对象，没有包可以直接用它
#arm相关操作
arms <- GenomicRanges::setdiff(chromosomes, tcmeres)
arms <- arms[-c(25,27,29,41,43)]##跳过

armlevels <- c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q",
               "7p","7q","8p","8q", "9p", "9q","10p","10q","11p","11q","12p",
               "12q","13q","14q","15q","16p","16q","17p","17q","18p","18q",
               "19p", "19q","20p","20q","21q","22q")##大致就是如何挑选全基因组有用的位点，我们用不到

arms$arm <- armlevels
AB <- AB[-queryHits(findOverlaps(AB, gaps.hg19))]
AB <- AB[queryHits(findOverlaps(AB, arms))]
AB$arm <- armlevels[subjectHits(findOverlaps(AB, arms))]##经过前面的操作，将AB这个表格增加了一列arm列

seqinfo(AB) <- seqinfo(Hsapiens)[seqlevels(seqinfo(AB))]#seqlevels得到的是每一个水平的名字，例如chr1这样的信息。
#相当于从Hsapien中取出AB中对应的信息，然后给AB添加信息,目的是让AB的信息更全。


#seqnames seqlengths isCircular genome
#  chr1      249250621      FALSE   hg19

AB <- trim(AB)#修剪越界的非NA的非环状序列
AB$gc <- GCcontent(Hsapiens, AB)#GCcontent(Hsapiens, GRanges("chr1", IRanges(1e6, 1e6 + 1000)))在某一个位置的区间GC含量

## These bins had no coverage
AB <- AB[-c(8780, 13665)]#-c(A,B)的是小于0的,这里选取了没reads覆盖的地方。
fragments <- readRDS(fragfile)#导入片段信息
# 
### Filters
fragments <- fragments[-queryHits(findOverlaps(fragments, filters.hg19))]#过滤掉之前已经从基因组过滤的区间片段,其中负号在R的index里面算作是排除符号
w.all <- width(fragments)#width指的是区间的长度

fragments <- fragments[which(w.all >= 100 & w.all <= 220)]#筛选符合片段长度要求的区间
w <- width(fragments)#得到fragment的每一个片段的长度,这里发现fragment是每一条序列的片段位置信息，一行指的是一条序列

frag.list <- split(fragments, w)#将每一种长度的区间进行区分,形成一个列表，列表每一个名字叫片段的长度

counts <- sapply(frag.list, function(x) countOverlaps(AB, x))#这里的overlap更多的指的是包含关系，直接涵盖的关系才能算作一个overlap,这里得到的是每一个区域overlap的个数
                 #指的是只需要在特定区间内的序列片段,最终返回的是每一个区间的个数
if(min(w) > 100) {
    m0 <- matrix(0, ncol=min(w) - 100, nrow=nrow(counts),#行为每一种片段长度，列为每一个区域的片段个数
                 dimnames=list(rownames(counts), 100:(min(w)-1)))
    counts <- cbind(m0, counts)#最后得到的是区间为行，前几列是100-最短长度，其余是最短到220的片段个数
}

olaps <- findOverlaps(fragments, AB)#得到序列和100kb区间的overlap位置
bin.list <- split(fragments[queryHits(olaps)], subjectHits(olaps))#得到列表，名称是fragment的序列号1，2，3，内容为AB的被overlap的片段信息
bingc <- rep(NA, length(bin.list))#多少个被匹配到的fragment就有多少个NA
bingc[unique(subjectHits(olaps))] <- sapply(bin.list, function(x) mean(x$gc))#对每一个fragment匹配上的位置求GC含量的均值，最后得到每一个fragment区段的gc含量均值

### Get modes
Mode <- function(x) {#该函数可以给出片段长度出现频率最高的一个片段区间
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]#tabulate为了看在固定区间内，每一个数出现了几次
}
modes <- Mode(w)##评选出出现频率最高的片段区间
medians <- median(w)#取出fragment长度的均值
q25 <- quantile(w, 0.25)#取出1/4处的长度
q75 <- quantile(w, 0.75)#取出3/4处的长度

short <- rowSums(counts[,1:51])#这个区间就是他们的到的结果100-150的
long <- rowSums(counts[,52:121])#这个区间是150-220的个数
ratio <- short/long#这个是区间每一个bin区间的ratio
short.corrected=gc.correct(short, bingc)#所有短片段进行GC覆盖度矫正
long.corrected=gc.correct(long, bingc)
nfrags.corrected=gc.correct(short+long, bingc)#对所有区间内的进行GC矫正
ratio.corrected=gc.correct(ratio, bingc)#对ratio进行GC矫正

AB$short <- short
AB$long <- long
AB$ratio <- short/long
AB$nfrags <- short+long
AB$short.corrected <- short.corrected
AB$long.corrected <- long.corrected
AB$nfrags.corrected <- nfrags.corrected
AB$ratio.corrected <- ratio.corrected

AB$mode <- modes
AB$mean <- round(mean(w), 2)
AB$median <- medians
AB$quantile.25 <- q25
AB$quantile.75 <- q75
AB$frag.gc <- bingc

for(i in 1:ncol(counts)) elementMetadata(AB)[,colnames(counts)[i]] <- counts[,i]

saveRDS(AB, filename)
q('no')
