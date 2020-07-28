library(tidyverse)
library(GenomicRanges)
bindir <- "../bins_100kb"

metadata <- read_csv("sample_reference.csv")
ids <- metadata %>% select(`WGS ID`) %>% unlist()#选取WGSid
files <- file.path(bindir, paste0(ids, "_bin_100kb.rds"))

bins.list <- lapply(files, readRDS)#打开每一个文件,并且将每一个文件打开之后的结果保存在一个list里面，成为bins.list
tib.list <- lapply(bins.list, as_tibble)##将文件类型改变为数据库table文件
names(tib.list) <- ids#给每一个列表的结果添加名字
#$this_name
# A tibble: 8,720 x 33
tib.list <- map2(tib.list, names(tib.list), ~ mutate(.x, id = .y)) %>%#mutate对dataframe的列名进行更改,这里给每一个tib文件创建了id列，该列为tib.list的names
    bind_rows() %>% select(id, everything())#bind_rows表示两个表或者向量上下拼接，在选取特定的列。

tib.list <- tib.list %>% select(-matches("X"))#
saveRDS(tib.list, "bins_100kbcompartments.rds")
