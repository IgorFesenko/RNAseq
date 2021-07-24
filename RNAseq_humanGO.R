# 4.3 GO ######
# устанавливем нужные пакеты
BiocManager::install("topGO")
BiocManager::install("Rgraphviz")
BiocManager::install("org.Hs.eg.db") # тут GO аннотация генома человека

library(topGO)

rm(list=ls())

qv = readRDS('/Users/igorfesenko/NGS/RNAseq/Rdata/ge.qv.Rds')
colnames(qv) [1] <- 'Qval'

s = factor(as.integer(apply(qv<0.05,1,sum)>0)) # вектор указывающий какие гены значимы (там стоит 1)
table(s) # 255 значимых гена
names(s) = rownames(qv) # имена указывают какие это гены. Незначимые гены (те у которых 0) составляют бэкграундный набор генов в сравнении с которым будет исследоваться обогащение значимых генов функциями. От этого набора сильно зависит результат. Он точно должен быть ограничен генами которые прошли фильтрацию и участвовали в тесте

# готовим данные для тесте
tgo1 = new("topGOdata", ontology = "BP", # будем смотреть на биологические функции, еще есть MF - молекулярная функция и CC - клеточный компонет
           allGenes = s,
           nodeSize = 5, # смотрим только на категории в которых есть хотябы 5 тестируемых генов
           annotationFun = annFUN.org,mapping='org.Hs.eg.db',ID='Ensembl') 

r1 = runTest(tgo1, algorithm = "classic", statistic = "fisher") # прогоняем тест
hist(score(r1)) # распределение p-value 


goqv = p.adjust(score(r1),m='BH') # коррекция на можественное тестирование (на количество категорий)
sort(goqv)[1:10] 
GenTable(tgo1,r1,topNodes=10) #10 самых значимых, во-всяком случае что-то нервное)

dev.off()
showSigOfNodes(tgo1, score(r1), firstSigNodes = 5, useInfo ='all') # нарисуем кусок графа GO со значимыми узлами
# также можно смотреть на другие группы генов - KEGG, reactome, EC, наличие эволюционных доменов, белковых фич (сайты модификаций, мембранные домены, неструктурированые регионы и т.д.)

