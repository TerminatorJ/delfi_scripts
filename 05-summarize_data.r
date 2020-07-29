library(tidyverse)
library(readxl)
#统计所有的样本，也可以用来构建基线
master <- read_csv("sample_reference.csv")
df.fr3 <- readRDS("../inst/extdata/bins_5mbcompartments.rds")
##这个是按照bin分类的结果
healthy.median <- df.fr3 %>%##构建基线，作为为正常人的样本，将每一个bin的中位数拿来做基线，因此median就是用来做基线的值，每一个bin分别对不同的样本进行统计
    group_by(bin) %>% ##对每一个5m bin 进行操作
    summarize(median.cov=median(nfrags2, na.rm=TRUE),##按照样本统计每一个指标在所有5m下的均值
              median.short=median(short2, na.rm=TRUE),
              median.long=median(long2, na.rm=TRUE),
              median.ratio=median(ratio2, na.rm=TRUE),
              median.corrected.cov=median(nfrags.corrected2, na.rm=TRUE),
              median.corrected.short=median(short.corrected2, na.rm=TRUE),
              median.corrected.long=median(long.corrected2, na.rm=TRUE),
              median.corrected.ratio=median(ratio.corrected2, na.rm=TRUE),#注意这里取的是中位数，而之前100kbbin里面取的是mean
              median.corrected.ratio2=median(short.corrected2/long.corrected2, na.rm=TRUE))#这里也将中位数的5mbin ration进行了一个
##这个是按照样本分类的结果，后缀为cor都是相关性的结果
summary.df <- df.fr3 %>% ungroup() %>% group_by(sample, type) %>%#对样本和疾病类型进行归类,即对health，lung， lymphocyte进行分别计算
    summarize(cov.cor=cor(nfrags2, healthy.median$median.cov, method="pearson", use="complete.obs"),#计算每一个样本与基线的每一个指标的pearson相关性
              short.cor=cor(short2, healthy.median$median.short, method="pearson", use="complete.obs"),
              long.cor=cor(long2, healthy.median$median.long, method="pearson", use="complete.obs"),
              ratio.cor=cor(ratio2, healthy.median$median.ratio, method="pearson", use="complete.obs"),
              cov.corrected.cor=cor(nfrags.corrected2, healthy.median$median.corrected.cov, method="pearson", use="complete.obs"),
              short.corrected.cor=cor(short.corrected2, healthy.median$median.corrected.short, method="pearson", use="complete.obs"),
              long.corrected.cor=cor(long.corrected2, healthy.median$median.corrected.long, method="pearson", use="complete.obs"),
              ratio.corrected.cor=cor(ratio.corrected2, healthy.median$median.corrected.ratio, method="pearson", use="complete.obs"),
              ratio2.corrected.cor=cor(short.corrected2/long.corrected2, healthy.median$median.corrected.ratio2, method="pearson", use="complete.obs"),
              nfrags = sum(nfrags2),
              mode_size=unique(mode_size),
              mean_size=unique(mean_size),
              median_size=unique(median_size),
              q25_size=unique(q25_size),
              q75_size=unique(q75_size),
              hqbases_analyzed = 100*sum(nfrags)*2,
              coverage = hqbases_analyzed/(504*5e6)
              )

summary.df <- inner_join(summary.df, master, by=c("sample"="WGS ID"))#将计算之后的值与样本信息结合起来
summary.df$`type` = relevel(as.factor(summary.df$`type`), "Healthy")#将所有的level都改成health

saveRDS(summary.df, "../inst/extdata/summary_tibble.rds")
