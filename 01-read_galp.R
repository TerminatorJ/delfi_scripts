library(GenomicAlignments)
library(GenomicRanges)
library(getopt)

### Used for getting information from shell script
args <- commandArgs(trailingOnly = TRUE)##仅仅返回--之后的字符串:[1] "--a"   "file1" "--b"   "file2" "--c"   "file3"
hh <- paste(unlist(args), collapse = " ")#将参数全部使用空格隔开形成一个大字符串:[1] "--a file1 --b file2 --c file3"
listoptions <- unlist(strsplit(hh, "--"))[-1]##去掉--：[1] "a file1 " "b file2 " "c file3" 
#得到每一个参数传入的变量
#a file1  b file2   c file3 
# "file1"  "file2"  "file3" 
options.args <- sapply(listoptions, function(x) {
    unlist(strsplit(x, " "))[-1]
})
##得到每一个参数名称
#a file1  b file2   c file3 
#     "a"      "b"      "c"
options.names <- sapply(listoptions, function(x) {
    option <- unlist(strsplit(x, " "))[1]
})
#得到每一个传入的参数以及其对应的名称
# a       b       c 
#"file1" "file2" "file3" 
names(options.args) <- unlist(options.names)
id <- options.args[1]
bamdir <- options.args[2]
galpdir <- options.args[3]
###

### Read GAlignmentPairs
bamfile <- file.path(bamdir, id)##得到路径类型的字符串a/b$
indexed.bam <- gsub("$", ".bai", bamfile)#给文件添加末尾添加.bai，$表示字符串的末尾,结果应该是bam.bai
if (!file.exists(indexed.bam)) {
    indexBam(bamfile)#对bam文件快速构建索引，产生bai文件，这之后，所有的bam文件都被加上索引
}
#ScanBamParam可以用来进行Bam文件的过滤，scanBamFlag可以用来指定特定的筛选条件
param <- ScanBamParam(flag = scanBamFlag(isDuplicate = FALSE,#去重
                                         isSecondaryAlignment = FALSE,#去除多比对序列
                                         isUnmappedQuery = FALSE),#去除没有比对上的序列
                      mapqFilter = 30)#不符合Q30的序列被去除掉
sample <- gsub(".bam", "", id)#id指的是每一个样本比对的文件，可以看出原来的命名为sample.bam

galp.file <- file.path(galpdir, paste0(sample, ".rds"))#设定保存文件的文件名称和路径
galp <- readGAlignmentPairs(bamfile, param = param)#将过滤后的文件进行读取
saveRDS(galp, galp.file)#存储过滤后的文件
