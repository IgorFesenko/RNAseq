BiocManager::install("edgeR") 
library(edgeR)  # загрузка библиотеки

rm(list=ls())

# загружаем данные которые сохранили в прошлый раз
counts = readRDS('/Users/igorfesenko/NGS/RNAseq/Rdata/filtered_counts.Rdata')
meta = readRDS('/Users/igorfesenko/NGS/RNAseq/Rdata/description.Rdata')

# Нормализуем данные
edger = DGEList(counts) # создаем объект edgeR хранящий каунты
edger = calcNormFactors(edger,method='RLE') # нормируем методом RLE


edger$samples # нормировочный факторы близки к единице, все в порядке
edger$dispersion # смотрим на дисперсию

#создаем матрицу дизайна. мы будем искать отличия между регионами мозга (tissue),
design = model.matrix(~ tissue,data = meta)
design

# определяем биологическую вариабельность
edger = estimateDisp(edger,design)

# рисуем зависимость биологической вариабельности от средней экспрессии, 
#немного падает с ростом экпсрессии
plotBCV(edger)

# строим модель
glm = glmFit(edger,design)

# считыем p-value для ткани
tissue = glmLRT(glm,2:4)
names(tissue)

# размер эффекта, средний уровень экспрессии, отношение правдоподобий и p-value
tissue$table[1:2,]

hist(tissue$table$PValue) # смотрим на распределения p-value

# сделаем поправку на множественное тестирование
pv = cbind(tissue=glmLRT(glm,2:4)$table$PValue)
rownames(pv) = rownames(counts)
qv = apply(pv,2,p.adjust,m='BH')
apply(qv < 0.05,2,sum) # количество значимо меняющихся генов


saveRDS(qv,'/Users/igorfesenko/NGS/RNAseq/Rdata/ge.qv.Rds') # сохраним корректированные p-value для последующего использования

# Выводим список 10 наиболее меняющихся генов
cpm = cpm(edger) # посчитаем cpm с учетом RLE нормировки
qv[order(qv[,1])[1:10],] # 10 самых значимых генов


############################
# 4.2 кластеризация #######
cpm.s = cpm[apply(qv<0.05,1,sum)>0,] # отбираем cpm для генов значимых хотябы в одном сравнении
cpm.s = t(scale(t(cpm.s))) # z-score
head(cpm.s)

hcl = hclust(as.dist(1-cor(t(cpm.s)))) # иерархическая кластеризация

plot(hcl)
cl = cutree(hcl,3) # режем дерево на 3 кластера
table(cl) # количетсво генов в кластерах

# кластеры приянто упорядочивать по количеству генов
# напишем функцию для этого
renameClustsBySize = function(x){
  t = table(x)
  n = names(t)
  o = order(t,decreasing=T)
  r = x
  for(i in 1:length(o))
    r[x == n[o[i]]] = i
  r
}
cl = renameClustsBySize(cl)
table(cl) #теперь упорядочины



