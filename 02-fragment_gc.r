tx <-as.numeric(Sys.getenv("SGE_TASK_ID"))#读取环境变量，环境变量可以在配置文件里面设置，假定为 ABCDEFG
galpdir <- "/dcl01/scharpf1/data/galignmentpairs/Low_Coverage_WGS"#假定是这个路径"/Users/wangjun/DELFI"

files <- list.files(galpdir, full.names=TRUE)#列出该路径下的所有文件,包含该文件的文件路径+文件名称
#[1] "/Users/wangjun/DELFI/AB_chr1.rds"                        
#[2] "/Users/wangjun/DELFI/chr1_tibble.rds"                    
#[3] "/Users/wangjun/DELFI/facer_grid.pdf" 
file <- files[tx]#得到目标文件。如：[1] "/Users/wangjun/DELFI/AB_chr1.rds"

names <- strsplit(basename(file), "_")[[1]][1]#basename获取路径末尾的文件名称,最终是获取的样本的名字
#[[1]]
#[1] "AB"       "chr1.rds"

outdir <- "/dcl01/scharpf1/data/galignmentpairs/Low_Coverage_WGS/fragments"
# if(file.exists(file.path(outdir, paste0(names, "_frags.rds")))) q('no')

library(GenomicAlignments)
library(GenomicRanges)
library(Rsamtools)
library(devtools)
library(Homo.sapiens)
library(biovizBase)
library(BSgenome.Hsapiens.UCSC.hg19)
class(Homo.sapiens)


galp <- readRDS(file)#导入目标文件
frags <- granges(keepSeqlevels(galp, paste0("chr", 1:22), pruning.mode="coarse"),#keepSeqlevels仅仅保留哪些seqnames，按照levels过滤，
             on.discordant.seqnames="drop")#之后创建Granges对象

## filter outliers
w.all <- width(frags)#得到每一条cfDNA的长度
q.all <- quantile(w.all, c(0.001, 0.999))#得到0.1%长度的cfDNA和99.9%的长度的cfDNA，这个区间以外都被认为是outliers
frags <- frags[which(w.all > q.all[1] & w.all < q.all[2])]#过滤区间之内的值

gcs <- GCcontent(Hsapiens, unstrand(frags))#获取每一个cfDNA片段的GC含量
frags$gc <- gcs#将GC含量添加到frags的其中一个列中

saveRDS(frags, file.path(outdir, paste0(names, "_frags.rds")) )
q('no')
