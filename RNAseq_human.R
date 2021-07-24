#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

#BiocManager::install("Rsubread")
rm(list=ls())

#загружаем библиотеки
library(Rsubread)

# читаем предобработанные данные из featureCounts
fc <- read.table('/Users/igorfesenko/read_counts.tsv', sep='\t', header=TRUE)
colnames(fc)
row.names(fc) <- fc$Geneid
fc$Geneid <- NULL

dim(fc) # размер таблицы гены*образцы
head(fc)


#загружаем таблицу с описанием данных

samples <- read.table('/Users/igorfesenko/samples_brain_R.csv', sep='\t', header=TRUE)
rownames(samples) = samples$sample
# выберем цвета которыми будем отрисовывать образцы - по ткани
samples$col = c('Prefrontal (BA10)'='orange',
             'Hypothalamus'='red',
             'Cerebellar White Matter'='green',
             'Cerebellar Grey Matter'='blue')[samples$tissue]
samples



# ФИЛЬТРАЦИЯ ГЕНОВ ПО ПОКРЫТИЮ

#распределение генов по суммарному окрытию ридами.
hist(log10(1+apply(fc,1,mean)),20)
abline(v=log10(1+10),col='red') # порог в 10 ридов в среднем на образец

# гены которыe вообще не экспрессируются в этих данных
table(apply(fc,1,mean) == 0)

# гены которые имеют покрытие больше порога - больше 10 ридов в среднем на образец
table(apply(fc,1,mean) >= 10) 

# возьмем только гены с покрытием больше 10 ридов на образец
fc10 = fc[apply(fc,1,mean) >= 10,]
table(apply(fc10,1,mean) >= 10)

# Отобрали 20890 генов для дальнейшего анализа
# сохраним обработанные данные на будущее
saveRDS(samples,'/Users/igorfesenko/NGS/RNAseq/Rdata/description.Rdata')
saveRDS(fc10,'/Users/igorfesenko/NGS/RNAseq/Rdata/filtered_counts.Rdata')

# 2.2 посмотрим немного на данные #######
# распределение не похоже на лог-нормальное
hist(fc10[,1])
hist(log(fc10[,1]))

# самосогласованность
cors = cor(fc10,m='sp') # посчитаем матрицу коэффициентов корреляции Спирмана между образцами (колонками)
corp = cor(fc10,m='p') # Пирсона
corlp = cor(log(1+fc10),m='p') # Пирсона в логшкале
cors # напечатаем матрицу

# построим хитмап
heatmap(1-cors, symm = T, distfun = function(x){as.dist(x)})

# MDS - построим для разных типов корреляций

mdss =  cmdscale(1-cors ,k = 2) # Спирмен, количество размерностей - 2
mdsp =  cmdscale(1-corp ,k = 2) # Пирсон
mdslp = cmdscale(1-corlp,k = 2) # Пирсон в логшкале

par(mfrow=c(2,2)) # разделим экран на 4 графика, чтобы нарисовать все MDS разом

plot(mdss ,pch=19,col=samples$col,cex=2,xlab='Dim 1',ylab='Dim 2',
     main='Spearman', bty='n')  
plot(mdsp, pch=19,col=samples$col,cex=2,xlab='Dim 1',ylab='Dim 2',main='Pearson'    ,bty='n')
plot(mdslp, pch=19,col=samples$col,cex=2,xlab='Dim 1',ylab='Dim 2',main='log Pearson',bty='n')


# PCA
#сначала нормируем на размер библиотеки
cpm = sweep(fc10,2,apply(fc10, 2,sum),'/')*1e6 # номируем данные на размер библиотеки.

pca = prcomp(t(cpm),scale. = TRUE)
barplot(pca$sdev) 
# количество вариебельности объясненное каждой последующей компонентой

# вычисляем долю объясненный вариабельности
var = paste0('PC',1:length(pca$sdev),' ',round(pca$sdev/sum(pca$sdev)*100,1),'%') 
plot(pca$x,xlab=var[1],ylab=var[2],bty='n',col=samples$col,cex=2,pch=19)
