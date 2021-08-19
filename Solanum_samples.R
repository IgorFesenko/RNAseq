# Загружаем данные по каунтам после featureCounts
rm(list=ls())

fc <- read.table('/Users/igorfesenko/Solanum_phytozome_read_counts', sep='\t', header=TRUE)
colnames(fc)
dim(fc)

# имена колонок содержать ".bam", уберем:
colnames(fc) <-  sub('.bam','',colnames(fc))
colnames(fc) <-  sub('bam.','',colnames(fc))

#сохраняем только данные по каунтам
filt_fc <-  fc[,c("Geneid", "rnaseq10_s","rnaseq11_s","rnaseq12_s","rnaseq13_s","rnaseq14_s",
                "rnaseq15_s","rnaseq16_s","rnaseq1_s","rnaseq2_s","rnaseq3_s","rnaseq4_s","rnaseq5_s",
                "rnaseq6_s","rnaseq7_s","rnaseq8_s","rnaseq9_s")]

#меняем порядок столбцов
colorder <-  c("Geneid","rnaseq1_s","rnaseq2_s","rnaseq3_s","rnaseq4_s","rnaseq5_s",
                       "rnaseq6_s","rnaseq7_s","rnaseq8_s","rnaseq9_s","rnaseq10_s","rnaseq11_s","rnaseq12_s","rnaseq13_s","rnaseq14_s",
                       "rnaseq15_s","rnaseq16_s")

filt_fc_col <- filt_fc[, colorder]

row.names(filt_fc_col) <- filt_fc_col$Geneid
filt_fc_col$Geneid <- NULL

head(filt_fc_col)
# сохраняем переменную counts в файл
saveRDS(filt_fc_col,'counts.Rdata')

# создадим табличку с информацией об образцах. В данном случае вся информация есть в именах образцов, иногда она находится в отдельных таблицах
meta <-  data.frame(genotype=c("mock","mock","mock","mock","mock","mock","mock","mock","PVY","PVY","PVY","PVY","PVY","PVY","PVY","PVY"),
                   temp=c(22,22,22,22,28,28,28,28,22,22,22,22,28,28,28,28),
                   col=c("green","green","green","green","orange","orange","orange","orange","blue","blue","blue","blue","red","red","red","red"))

rownames(meta) = colnames(filt_fc_col)

head(meta)

# ФИЛЬТРАЦИЯ ГЕНОВ ПО ПОКРЫТИЮ

#распределение генов по суммарному окрытию ридами.
hist(log10(1+apply(filt_fc_col,1,mean)),20)
abline(v=log10(1+10),col='red') # порог в 10 ридов в среднем на образец

# гены которыe вообще не экспрессируются в этих данных
table(apply(filt_fc_col,1,mean) == 0)

# гены которые имеют покрытие больше порога - больше 10 ридов в среднем на образец
table(apply(filt_fc_col,1,mean) >= 10) 

# возьмем только гены с покрытием больше 10 ридов на образец
fc10 = filt_fc_col[apply(filt_fc_col,1,mean) >= 10,]
table(apply(fc10,1,mean) >= 10)

# Отобрали 24517 генов для дальнейшего анализа
# сохраним обработанные данные на будущее
saveRDS(fc10,'/Users/igorfesenko/Google Диск/solanum/RNAseq/Rdata/filtered_counts.Rdata')
saveRDS(meta,'/Users/igorfesenko/Google Диск/solanum/RNAseq/Rdata/meta.Rdata')


#EDA
hist(fc10[,1])
hist(log(fc10[,3]))

plot(filt_fc_col[,1],filt_fc_col[,2],log='xy')

#сравним все образцы
pairs(filt_fc_col,log='xy',pch='.')

# самосогласованность
cors = cor(fc10,m='sp') # посчитаем матрицу коэффициентов корреляции Спирмана между образцами (колонками)
corp = cor(fc10,m='p') # Пирсона
corlp = cor(log(1+fc10),m='p') # Пирсона в логшкале
cors # напечатаем матрицу

# построим хитмап
heatmap(1-cors, symm = T, margins=c(5,5),distfun = function(x){as.dist(x)})

# MDS - построим для разных типов корреляций

mdss =  cmdscale(1-cors ,k = 2) # Спирмен, количество размерностей - 2
mdsp =  cmdscale(1-corp ,k = 2) # Пирсон
mdslp = cmdscale(1-corlp,k = 2) # Пирсон в логшкале

par(mfrow=c(2,2)) # разделим экран на 4 графика, чтобы нарисовать все MDS разом

plot(mdss ,pch=19,col=meta$col,cex=2,xlab='Dim 1',ylab='Dim 2',
     main='Spearman', bty='n')
plot(mdsp, pch=19,col=meta$col,cex=2,xlab='Dim 1',ylab='Dim 2',main='Pearson'    ,bty='n')
plot(mdslp, pch=19,col=meta$col,cex=2,xlab='Dim 1',ylab='Dim 2',main='log Pearson',bty='n')

# PCA
#сначала нормируем на размер библиотеки
cpm = sweep(fc10,2,apply(fc10, 2,sum),'/')*1e6 # номируем данные на размер библиотеки.

pca = prcomp(t(cpm),scale. = TRUE)
barplot(pca$sdev) 
# количество вариебельности объясненное каждой последующей компонентой

# вычисляем долю объясненный вариабельности
var = paste0('PC',1:length(pca$sdev),' ',round(pca$sdev/sum(pca$sdev)*100,1),'%') 
plot(pca$x,xlab=var[1],ylab=var[2],bty='n',col=meta$col,cex=2,pch=19)
