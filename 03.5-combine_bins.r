#将所有的样本100kb的bin结合起来形成一个大表格
library(tidyverse)
library(GenomicRanges)
bindir <- "../bins_100kb"
##目的，产生一个list，每一个键包含WGS_id的一个号码，其实相当于每一个样本，每一个样本得到一个AB列表，之后就可以对其进行批量处理

metadata <- read_csv("sample_reference.csv")
ids <- metadata %>% select(`WGS ID`) %>% unlist()#选取WGSid
files <- file.path(bindir, paste0(ids, "_bin_100kb.rds"))

bins.list <- lapply(files, readRDS)#打开每一个文件,并且将每一个文件打开之后的结果保存在一个list里面，成为bins.list
tib.list <- lapply(bins.list, as_tibble)##将文件类型改变为数据库table文件
names(tib.list) <- ids#给每一个列表的结果添加名字
#$this_name
# A tibble: 8,720 x 33
#给tib.list的每一个列表里面的内容创建了一个列，为其WGS的id名
tib.list <- map2(tib.list, names(tib.list), ~ mutate(.x, id = .y)) %>%#mutate对dataframe的列名进行更改,这里给每一个tib文件创建了id列，该列为tib.list的names
    bind_rows() %>% select(id, everything())#bind_rows表示两个表或者向量上下拼接，在选取特定的列。select(id,everything())表示将id放在最前头
#样例：
#$this_name#有很多的键，每一个都执行以下操作
## A tibble: 8,720 x 33
#   sample type  chromosome arm   combine short2 long2 short.corrected2 long.corrected2 hic.eigen id
#   <chr>  <fct> <fct>      <fct>   <dbl>  <dbl> <dbl>            <dbl>           <dbl>     <dbl>  1
# 1 PGDX1… Lymp… chr1       1p          1  23289 29487           14419.          29010.    -1.61   1
#最后形成的是一个大表格，包含了所有样本的bin信息，同事添加了id列

tib.list <- tib.list %>% select(-matches("X"))#将X列排除掉
saveRDS(tib.list, "bins_100kbcompartments.rds")
