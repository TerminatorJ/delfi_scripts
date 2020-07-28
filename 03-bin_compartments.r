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
load("./filters.hg19.rda")#被分为1565个bin，该文件是一个GRanges对象,需要明确这些bin是如何被过滤的，这个表示profile图片上的横坐标

library(RCurl)##这个包用来访问外部网络，进行下载东西使用，到那时这里我们提前下载到本地了。
ABurl <- getURL('https://raw.githubusercontent.com/Jfortin1/HiC_AB_Compartments/master/data/hic_compartments_100kb_ebv_2014.txt', ssl.verifyhost=FALSE, ssl.verifypeer=FALSE)
#因此我们直接导入本地包：
AB=read.table("hic_compartments_100kb_ebv_2014.txt",header=TRUE)##从文件中可以看到，作者将区域的bin划分为100kb，在图上显示的是一个点，所有的点结合起来就是折线
AB <- read.table(textConnection(ABurl), header=TRUE)
AB <- makeGRangesFromDataFrame(AB, keep.extra.columns=TRUE)#将其变成Granges对象

chromosomes <- GRanges(paste0("chr", 1:22),
                       IRanges(0, seqlengths(Hsapiens)[1:22]))#构建Granges对象，包含的是hg19的每一个信息
#seqnames      ranges strand
#          <Rle>   <IRanges>  <Rle>
#   [1]     chr1 0-249250621      *
#   [2]     chr2 0-243199373      *
#   [3]     chr3 0-198022430      *
#   [4]     chr4 0-191154276      *
#   [5]     chr5 0-180915260      *
#   ...      ...         ...    ...
#  [18]    chr18  0-78077248      *
#  [19]    chr19  0-59128983      *
#  [20]    chr20  0-63025520      *
#  [21]    chr21  0-48129895      *
#  [22]    chr22  0-51304566      *
#  -------

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
#相当于从Hsapien中取出AB中对应的信息，然后给AB添加信息,目的是让AB的信息更全,  经过比较发现，这样去比对其实没有意义，因为序列的长度在AB和Hsapiens中是一样的
#seqnames seqlengths isCircular genome
#  chr1      249250621      FALSE   hg19
#  chr2      243199373      FALSE   hg19
#  chr3      198022430      FALSE   hg19
#  chr4      191154276      FALSE   hg19
#  chr5      180915260      FALSE   hg19
#  ...             ...        ...    ...
#  chr18      78077248      FALSE   hg19
#  chr19      59128983      FALSE   hg19
#  chr20      63025520      FALSE   hg19
#  chr21      48129895      FALSE   hg19
#  chr22      51304566      FALSE   hg19

AB <- trim(AB)#修剪越界的非NA的非环状序列,确实存在这么多的修剪
#which(AB2!=AB)
# [1]  2204  2205  4551  6498  8369 11796 13326 19722 20680 22336 23092 23859 23860
#[14] 24600 25158
#必须要进过trim否则后面GCcount的时候会报错
#AB$gc <- GCcontent(Hsapiens, AB)
#Error in loadFUN(x, seqname, ranges) : 
#  trying to load regions beyond the boundaries of non-circular sequence "chr1"

AB$gc <- GCcontent(Hsapiens, AB)#GCcontent(Hsapiens, GRanges("chr1", IRanges(1e6, 1e6 + 1000)))在某一个位置的区间GC含量,得到了每一个100kbbin的区间gc含量

## These bins had no coverage
AB <- AB[-c(8780, 13665)]#这里排除了没reads覆盖的地方。
fragments <- readRDS(fragfile)#导入片段信息，目前没有找到片段信息，能确定其是一个GRanges对象
# 
### Filters
fragments <- fragments[-queryHits(findOverlaps(fragments, filters.hg19))]#过滤掉之前已经从基因组过滤的区间片段,其中负号在R的index里面算作是排除符号
w.all <- width(fragments)#width指的是区间的长度


fragments <- fragments[which(w.all >= 100 & w.all <= 220)]#筛选符合片段长度要求的区间
w <- width(fragments)#。统计每一个fragment的片段长度
#[1] 10  9  8  7  6  5  4  3  2  1
frag.list <- split(fragments, w)#将片段使用长度进行分裂，GRangesList列表的每一个名称为长度，内容为该长度下的所有片段

counts <- sapply(frag.list, function(x) countOverlaps(AB, x))#返回每一个bin里面有多少条序列
#     1 2 3 4 5 6 7 8 9 10
#    [1,] 0 0 0 0 0 0 0 0 0  0
#    [2,] 0 0 0 0 0 0 0 0 0  0
#    [3,] 0 0 0 0 0 0 0 0 0  0
#列名代表长度，行名代表每一个bin，得到每一个bin下面有多少条reads，并且知道这些reads的长度
if(min(w) > 100) {
    m0 <- matrix(0, ncol=min(w) - 100, nrow=nrow(counts),#行为每一种片段长度，列为每一个区域的片段个数
                 dimnames=list(rownames(counts), 100:(min(w)-1)))
    counts <- cbind(m0, counts)#最后得到的是区间为行，第一列从100开始，为了后面直接取相应的值方便，因为可能出现直接从120开始，但是后面取的时候取1：50（100-150）会出现越界的情况
}

olaps <- findOverlaps(fragments, AB)#得到cfDNA片段和bin匹配的位置
#      queryHits subjectHits       queryHits指的是fragments的位置；subjectHit指的是bin的位置
#      <integer>   <integer>
#  [1]         2        2205
#  [2]         3        2205
#  [3]         4        2205
#  [4]         7        4551
#  [5]         8        4551
#  [6]         9        4551
#  [7]        10        4551  
                 
bin.list <- split(fragments[queryHits(olaps)], subjectHits(olaps))#得到列表，得到每一个bin里面包含的所有fragment信息
#$`2205`
#GRanges object with 3 ranges and 2 metadata columns:
#    seqnames    ranges strand |     score        GC
#       <Rle> <IRanges>  <Rle> | <integer> <numeric>
#  b     chr2      2-10      + |         2  0.888889
#  c     chr2      3-10      + |         3  0.777778
#  d     chr2      4-10      * |         4  0.666667
#  -------
#  seqinfo: 3 sequences from an unspecified genome; no seqlengths

#$`4551`
#GRanges object with 4 ranges and 2 metadata columns:
#    seqnames    ranges strand |     score        GC
#       <Rle> <IRanges>  <Rle> | <integer> <numeric>
#  g     chr3      7-10      + |         7  0.333333
#  h     chr3      8-10      + |         8  0.222222
#  i     chr3      9-10      - |         9  0.111111
#  j     chr3        10      - |        10  0.000000
#  -------
#  seqinfo: 3 sequences from an unspecified genome; no seqlengths
                 
bingc <- rep(NA, length(bin.list))#多少个被匹配到的就有多少个NA,相当于创建新的数组
bingc[unique(subjectHits(olaps))] <- sapply(bin.list, function(x) mean(x$gc))#得到每一个bin的GC含量（均值）

### Get modes
Mode <- function(x) {#该函数可以给出片段长度出现频率最高的一个片段区间
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]#tabulate为了看在固定区间内，每一个数出现了几次
}
modes <- Mode(w)##评选出出现频率最高的片段区间
medians <- median(w)#取出fragment长度的均值
q25 <- quantile(w, 0.25)#取出fragment1/4处的长度
q75 <- quantile(w, 0.75)#取出fragment3/4处的长度

short <- rowSums(counts[,1:51])#这个区间就是他们的到的结果100-150的，前面已经进行过100的填充处理
long <- rowSums(counts[,52:121])#这个区间是150-220的个数
ratio <- short/long#这个是区间每一个bin区间的ratio
short.corrected=gc.correct(short, bingc)#所有短片段进行GC覆盖度矫正
long.corrected=gc.correct(long, bingc)#
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
