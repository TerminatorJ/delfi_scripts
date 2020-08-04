library(tidyverse)
library(multidplyr)
library(GenomicRanges)
library(readxl)

df.fr <- readRDS("bins_100kbcompartments.rds")#导入的是一个大table，table中添加了一列id来进行不同样本的区分
master <- read_csv("sample_reference.csv")#所有患者信息

df.fr2 <- inner_join(df.fr, master, by=c("id"="WGS ID"))##形成了一个大表格，包括样本信息也包括ratio信息

armlevels <- c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q",
               "7p","7q","8p","8q", "9p", "9q","10p","10q","11p","11q","12p",
               "12q","13q","14q","15q","16p","16q","17p","17q","18p","18q",
               "19p", "19q","20p","20q","21q","22q")
df.fr2$arm <- factor(df.fr2$arm, levels=armlevels)#得到了所有的需要的染色体臂，也包括过滤掉依稀而不用的臂

## combine adjacent 100kb bins to form 5mb bins. We count starting from
## the telomeric(端粒) end and remove the bin closest to the centromere（着丝粒） if it is
## smaller than 5mb.端粒开始计数，包含着丝粒的不足5mb的位点删除
#对df.fr添加一列conbine，将所有的染色体p进行分组，分成175段，其中每一段有50个值每一个值相当于100kb，100kb*50=5M的长度，这样就可以将样本分为5M的长度区间了
df.fr2 <- df.fr2 %>% group_by(id, arm) %>%
    mutate(combine = ifelse(grepl("p", arm), ceiling((1:length(arm))/50),#其实就是将某一个样本的8000多个bin使用ceiling将所有的arm分成175段,其中每一段有50个相同的值
                           ceiling(rev((1:length(arm))/50) )))#ceiling函数是向上取整函数，得到一个数值型向量的所有天花板整数。将函数得到的mutate添加一列在总表上面。
#对不同id，不同100kb的bin，不同arm，不同5m的bin进行分别计算
df.fr3 <- df.fr2 %>% group_by(id, seqnames, arm, combine) %>%
    summarize(short2=sum(short),
              long2=sum(long),
              short.corrected2=sum(short.corrected),
              long.corrected2=sum(long.corrected),
              hic.eigen=mean(eigen),
              gc=mean(C.G),
              ratio2=mean(ratio),
              ratio.corrected2=mean(ratio.corrected),##最终使用的还是5m的bin展示吗
              nfrags2=sum(nfrags),
              nfrags.corrected2=sum(nfrags.corrected),
              domain = median(as.integer(domain)),
              short.var=var(short.corrected),
              long.var=var(long.corrected),
              nfrags.var=var(nfrags.corrected),
              mode_size=unique(mode),
              mean_size=unique(mean),
              median_size=unique(median),
              q25_size=unique(quantile.25),
              q75_size=unique(quantile.75),
              start=start[1],
              end=rev(end)[1],
              binsize = n())
### assign bins
df.fr3 <- inner_join(df.fr3, master, by=c("sample"="WGS ID"))#得到所有样本的5m信息,由于summarise之后形成的表格没有样本信息，所以innerjoin一下来得到样本信息
df.fr3 <- df.fr3 %>% mutate(type = gsub(" Cancer|carcinoma", "", `Patient Type`, ignore.case=TRUE))#对每一个patient type进行修饰
df.fr3 <- df.fr3 %>% filter(binsize==50)#仅仅选取binsize在50个的bin舍弃掉着丝粒上面的bin
df.fr3 <- df.fr3 %>% group_by(sample) %>% mutate(bin = 1:length(sample))#对新形成的5m的bin进行排列，每一个bin代表一个5m的区间。

saveRDS(df.fr3, "bins_5mbcompartments.rds")
