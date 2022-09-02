# #PAQUETES NECESARIOS para todo el script (recomendable cargarlos al inicio)# #
library(ballgown)
library(NormalyzerDE)
library(FactoMineR)
library(factoextra)
library(limma)
library(clusterProfiler)
library(org.At.tair.db)
library(VennDiagram)
library("WGCNA") 
library("cluster")


# # # # # # # # # # # # # # # # # # # # # # # # # #
# EXPLORACIÓN INICIAL DE LOS DATOS Y NORMALIZACIÓN # 
# # # # # # # # # # # # # # # # # # # # # # # # # #

library(ballgown)
experimental.design <- read.csv("ecofun_experimental_design.csv")
experimental.design

ecofun.data <- ballgown(dataDir = ".", samplePattern = "sample", 
                        pData=experimental.design)

gene.expression <- gexpr(ecofun.data)
head(gene.expression)

table(experimental.design$ecotype)

colnames(gene.expression) <- c(
  "ex1velm_1",
  "ex1velm_2",
  "ex1velm_3",
  "ex1velt_1",
  "ex1velt_2",
  "ex1velt_3",
  "ex2velm_1",
  "ex2velm_2",
  "ex2velm_3",
  "ex2velm_4",
  "ex2velt_1",
  "ex2velt_2",
  "ex2velt_3",
  "ex3velm_1",
  "ex3velm_2",
  "ex3velt_1",
  "ex3velt_2",
  "ex3velt_3",
  "ex4velm_1",
  "ex4velm_2",
  "ex4velm_3",
  "ex4velt_1",
  "ex4velt_2",
  "ex4velt_3",
  "ex1bonm_1",
  "ex1bonm_2",
  "ex1bonm_3",
  "ex1bont_1",
  "ex1bont_2",
  "ex1bont_3",
  "ex2bonm_1",
  "ex2bonm_2",
  "ex2bonm_3",
  "ex2bont_1",
  "ex2bont_2",
  "ex2bont_3",
  "ex3bonm_1",
  "ex3bonm_2",
  "ex3bonm_3",
  "ex3bont_1",
  "ex3bont_2",
  "ex3bont_3",
  "ex1mojm_1",
  "ex1mojm_2",
  "ex1mojm_3",
  "ex1mojt_1",
  "ex1mojt_2",
  "ex1mojt_3",
  "ex2mojm_1",
  "ex2mojm_2",
  "ex2mojm_3",
  "ex2mojt_1",
  "ex2mojt_2",
  "ex2mojt_3",
  "ex3mojm_1",
  "ex3mojm_2",
  "ex3mojm_3",
  "ex3mojt_1",
  "ex3mojt_2",
  "ex3mojt_3",
  "ex1col0m_1",
  "ex1col0m_2",
  "ex1col0m_3",
  "ex1col0t_1",
  "ex1col0t_2",
  "ex1col0t_3",
  "ex2col0m_1",
  "ex2col0m_2",
  "ex2col0m_3",
  "ex2col0t_1",
  "ex2col0t_2",
  "ex2col0t_3",
  "ex3col0m_1",
  "ex3col0m_2",
  "ex3col0m_3",
  "ex3col0t_1",
  "ex3col0t_2",
  "ex3col0t_3",
  "ex4col0m_1",
  "ex4col0m_2",
  "ex4col0m_3",
  "ex4col0t_1",
  "ex4col0t_2",
  "ex4col0t_3")

write.table(x = gene.expression,file = "ecofun_gene_expression.tsv",
            quote = F,sep = "\t")

ecofun.gene.expression <- read.table(file = "ecofun_gene_expression.tsv",
                                     header = T,sep = "\t",as.is = T)

head(ecofun.gene.expression)
dim(ecofun.gene.expression)

png(filename = "raw_data_boxplot.png",width = 2000)
boxplot(log2(ecofun.gene.expression + 1),
        col=rainbow(ncol(ecofun.gene.expression)),ylab="log2(FPKM + 1)",
        cex.lab=1.5,las=2)
dev.off()

## Remove column 10 corresponding to ex3velm3 (degraded)
ecofun.gene.expression <- ecofun.gene.expression[,-10]

## Normalization. This step may be problematic as we are
## forcing all samples to have similar global distribution when
## major developmental processes are taking place
library(NormalyzerDE)

ecofun.gene.expression.df <- data.frame(rownames(ecofun.gene.expression),
                                        ecofun.gene.expression+1)
head(ecofun.gene.expression.df)
colnames(ecofun.gene.expression.df)[1] <- "geneID"
write.table(x = ecofun.gene.expression.df,file = "ecofun_gene_expression_df.tsv",
            quote = F,row.names = F,
            sep = "\t")

## Remove sample 10 from experimental design corresponding to ex3velm3
design <- data.frame(sample=colnames(ecofun.gene.expression),
                     group=experimental.design$ecotype[-10])

write.table(x = design,file = "ecofun_design.tsv",quote = F,row.names = F,
            sep = "\t")

normalyzer(jobName = "ecofun",designPath = "ecofun_design.tsv",
           dataPath = "ecofun_gene_expression_df.tsv",outputDir = ".")


#leer la matriz de expresión génica normalizada y transformada por log2. 
log.normalized.ecofun.gene.expression <- read.table(file="ecofun/Quantile-normalized.txt",
                                                    header=T,as.is=T)
head(log.normalized.ecofun.gene.expression)
geneID <- log.normalized.ecofun.gene.expression$geneID
dim(log.normalized.ecofun.gene.expression)
#Quédate sólo con las columnas desde la segunda a la 83 
log.normalized.ecofun.gene.expression <- as.matrix(log.normalized.ecofun.gene.expression[,2:84])
#nombra las filas con los nombres en la columna geneID.
rownames(log.normalized.ecofun.gene.expression) <- geneID
#Ya salta a la línea 329 y comienza el análisis de expresión diferencial.


#comprueba que los datos dentro del objeto son números 
is.numeric(log.normalized.ecofun.gene.expression)

#figura boxplot de datos normalizados
png(filename = "normalized_data_boxplot.png",width = 2000)
boxplot(log.normalized.ecofun.gene.expression,
        col=rainbow(ncol(log.normalized.ecofun.gene.expression)),
        ylab="log2(FPKM + 1)",cex.lab=1.5,las=2)
dev.off()


#ANÁLISIS DE COMPONENTES PRINCIPALES
library(FactoMineR)
library(factoextra)

pca.ecofun.gene.expression <- data.frame(colnames(log.normalized.ecofun.gene.expression),
                                         t(log.normalized.ecofun.gene.expression))
colnames(pca.ecofun.gene.expression)[1] <- "Ecotype"

res.ecofun.pca <- PCA(pca.ecofun.gene.expression, graph = FALSE,
                      scale.unit = TRUE,quali.sup = 1 )
res.ecofun.hcpc <- HCPC(res.ecofun.pca, graph=FALSE,nb.clust = -1)  
png(filename = "hcpc_ecofun.png",width = 2000,height = 1000)
fviz_dend(res.ecofun.hcpc,k=2,
          cex = 1,                       # Label size
          palette = "jco",               # Color palette see ?ggpubr::ggpar
          rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
          rect_border = "jco",           # Rectangle color
          type="rectangle",
          labels_track_height = 2400      # Augment the room for labels
)
dev.off()
png(filename = "pca_ecofun.png",width = 1000,height = 500)
fviz_pca_ind(res.ecofun.pca, col.ind = experimental.design$ecotype[-10], 
             pointsize=2, pointshape=21,fill="black",
             repel = TRUE, 
             addEllipses = TRUE,ellipse.type = "confidence",
             legend.title="Conditions",
             title="",
             show_legend=TRUE,show_guide=TRUE) 
dev.off()

gene.contribution.pca <- res.ecofun.pca$var$contrib
dim(gene.contribution.pca)

removeX <- function(x)
{
  return(as.numeric(substr(x=x,start = 2,stop = nchar(x))))  
}

sum(gene.contribution.pca[,"Dim.1"])

write(x = geneID[sapply(X = names(sort(gene.contribution.pca[,"Dim.1"],
                                       decreasing = T)[1:1000]),
                        FUN = removeX)],file = "1000_dim1.txt")

write(x = geneID[sapply(X = names(sort(gene.contribution.pca[,"Dim.2"],
                                       decreasing = T)[1:1000]),
                        FUN = removeX)],file = "1000_dim2.txt")


geneID[20669]

gene.coord.pca <- abs(res.ecofun.pca$var$coord)

write(x = geneID[sapply(X = names(sort(gene.coord.pca[,"Dim.1"],
                                       decreasing = T)[1:100]),
                        FUN = removeX)],file = "coord_100_dim1.txt")




normalized.ecofun.gene.expression <- (2^log.normalized.ecofun.gene.expression)-1

gene <- "AT2G23380"
gene.name <- "CLF"

ex1.col0.m <- normalized.ecofun.gene.expression[gene,paste0("ex1col0m_",1:3)]
ex1.col0.t <- normalized.ecofun.gene.expression[gene,paste0("ex1col0t_",1:3)]
ex1.vel.m <- normalized.ecofun.gene.expression[gene,paste0("ex1velm_",1:3)]
ex1.vel.t <- normalized.ecofun.gene.expression[gene,paste0("ex1velt_",1:3)]
ex1.moj.m <- normalized.ecofun.gene.expression[gene,paste0("ex1mojm_",1:3)]
ex1.moj.t <- normalized.ecofun.gene.expression[gene,paste0("ex1mojt_",1:3)]
ex1.bon.m <- normalized.ecofun.gene.expression[gene,paste0("ex1bonm_",1:3)]
ex1.bon.t <- normalized.ecofun.gene.expression[gene,paste0("ex1bont_",1:3)]

ex1.means <- c(mean(ex1.col0.m),mean(ex1.col0.t),
               mean(ex1.vel.m),mean(ex1.vel.t),
               mean(ex1.moj.m),mean(ex1.moj.t),
               mean(ex1.bon.m),mean(ex1.bon.t))
ex1.sds <- c(sd(ex1.col0.m),sd(ex1.col0.t),
             sd(ex1.vel.m),sd(ex1.vel.t),
             sd(ex1.moj.m),sd(ex1.moj.t),
             sd(ex1.bon.m),sd(ex1.bon.t))

ex2.col0.m <- normalized.ecofun.gene.expression[gene,paste0("ex2col0m_",1:3)]
ex2.col0.t <- normalized.ecofun.gene.expression[gene,paste0("ex2col0t_",1:3)]
ex2.vel.m <- normalized.ecofun.gene.expression[gene,paste0("ex2velm_",1:3)]
ex2.vel.t <- normalized.ecofun.gene.expression[gene,paste0("ex2velt_",1:3)]
ex2.moj.m <- normalized.ecofun.gene.expression[gene,paste0("ex2mojm_",1:3)]
ex2.moj.t <- normalized.ecofun.gene.expression[gene,paste0("ex2mojt_",1:3)]
ex2.bon.m <- normalized.ecofun.gene.expression[gene,paste0("ex2bonm_",1:3)]
ex2.bon.t <- normalized.ecofun.gene.expression[gene,paste0("ex2bont_",1:3)]

ex2.means <- c(mean(ex2.col0.m),mean(ex2.col0.t),
               mean(ex2.vel.m),mean(ex2.vel.t),
               mean(ex2.moj.m),mean(ex2.moj.t),
               mean(ex2.bon.m),mean(ex2.bon.t))

ex2.sds <- c(sd(ex2.col0.m),sd(ex2.col0.t),
             sd(ex2.vel.m),sd(ex2.vel.t),
             sd(ex2.moj.m),sd(ex2.moj.t),
             sd(ex2.bon.m),sd(ex2.bon.t))


ex3.col0.m <- normalized.ecofun.gene.expression[gene,paste0("ex3col0m_",1:3)]
ex3.col0.t <- normalized.ecofun.gene.expression[gene,paste0("ex3col0t_",1:3)]
ex3.vel.m <- normalized.ecofun.gene.expression[gene,paste0("ex3velm_",1:2)]
ex3.vel.t <- normalized.ecofun.gene.expression[gene,paste0("ex3velt_",1:3)]
ex3.moj.m <- normalized.ecofun.gene.expression[gene,paste0("ex3mojm_",1:3)]
ex3.moj.t <- normalized.ecofun.gene.expression[gene,paste0("ex3mojt_",1:3)]
ex3.bon.m <- normalized.ecofun.gene.expression[gene,paste0("ex3bonm_",1:3)]
ex3.bon.t <- normalized.ecofun.gene.expression[gene,paste0("ex3bont_",1:3)]

ex3.means <- c(mean(ex3.col0.m),mean(ex3.col0.t),
               mean(ex3.vel.m),mean(ex3.vel.t),
               mean(ex3.moj.m),mean(ex3.moj.t),
               mean(ex3.bon.m),mean(ex3.bon.t))

ex3.sds <- c(sd(ex3.col0.m),sd(ex3.col0.t),
             sd(ex3.vel.m),sd(ex3.vel.t),
             sd(ex3.moj.m),sd(ex3.moj.t),
             sd(ex3.bon.m),sd(ex3.bon.t))


ex4.col0.m <- normalized.ecofun.gene.expression[gene,paste0("ex4col0m_",1:3)]
ex4.col0.t <- normalized.ecofun.gene.expression[gene,paste0("ex4col0t_",1:3)]
ex4.vel.m <- normalized.ecofun.gene.expression[gene,paste0("ex4velm_",1:3)]
ex4.vel.t <- normalized.ecofun.gene.expression[gene,paste0("ex4velt_",1:3)]

ex4.means <- c(mean(ex4.col0.m),mean(ex4.col0.t),
               mean(ex4.vel.m),mean(ex4.vel.t),
               NA,NA,
               NA,NA)

ex4.sds <- c(sd(ex4.col0.m),sd(ex4.col0.t),
             sd(ex4.vel.m),sd(ex4.vel.t),
             NA,NA,
             NA,NA)


png(filename = paste(c(gene,"_",gene.name,".png"),collapse = ""),
    width = 1000)
barplot.data <- matrix(data = c(ex1.means,ex2.means,ex3.means,ex4.means),ncol=4)
colnames(barplot.data) <- paste0("Ex", 1:4)
colnames(barplot.data) <- c("30/12/20", "03/02/21", "18/02/21", "04/03/21")
rownames(barplot.data) <- c("Col0 Morning",
                            "Col0 Evening",
                            "Vel Morning",
                            "Vel Evening",
                            "Moj Morning",
                            "Moj Evening",
                            "Bon Morning",
                            "Bon Evening")


y.max <- max(c(ex1.means + ex1.sds,
               ex2.means + ex2.sds,
               ex3.means + ex3.sds,
               ex4.means + ex4.sds),na.rm = T)


xpos <- barplot(height = barplot.data,beside=T,col = c("red","lightpink1",
                                                       "blue","lightblue",
                                                       "darkgreen","lightgreen",
                                                       "brown","navajowhite3"),
                ylim=c(0,y.max),legend=rownames(barplot.data),
                main=paste(gene,gene.name,sep=" - "), cex.main=2,
                ylab="FPKM",cex.lab=2)

arrows(x0 = xpos[,1],y0 = ex1.means + ex1.sds,
       x1 = xpos[,1], y1 = ex1.means - ex1.sds,
       length = 0.02,angle = 90,code = 3)

arrows(x0 = xpos[,2],y0 = ex2.means + ex2.sds,
       x1 = xpos[,2], y1 = ex2.means - ex2.sds,
       length = 0.02,angle = 90,code = 3)

arrows(x0 = xpos[,3],y0 = ex3.means + ex3.sds,
       x1 = xpos[,3], y1 = ex3.means - ex3.sds,
       length = 0.02,angle = 90,code = 3)

arrows(x0 = xpos[,4],y0 = ex4.means + ex4.sds,
       x1 = xpos[,4], y1 = ex4.means - ex4.sds,
       length = 0.02,angle = 90,code = 3)
dev.off()


# # # # # # # # # # # # # # # # # # # # # 
# # ANÁLISIS DE EXPRESIÓN DIFERENCIAL # #
# # # # # # # # # # # # # # # # # # # # #

#leer la matriz de expresión génica normalizada y transformada por log2. 
log.normalized.ecofun.gene.expression <- read.table(file="ecofun/Quantile-normalized.txt",
                                                    header=T,as.is=T)
head(log.normalized.ecofun.gene.expression)
geneID <- log.normalized.ecofun.gene.expression$geneID
dim(log.normalized.ecofun.gene.expression)
#Quédate sólo con las columnas desde la segunda a la 83 
log.normalized.ecofun.gene.expression <- as.matrix(log.normalized.ecofun.gene.expression[,2:84])
#nombra las filas con los nombres en la columna geneID.
rownames(log.normalized.ecofun.gene.expression) <- geneID

#Es necesario caragar estos paquetes antes de correr la función
library(clusterProfiler)
library(org.At.tair.db)



##############################################################################
##############################################################################
################                                            ##################
################ FUNCIÓN EXTRAER CONTRASTES INDIVIDUALES    ##################
################                                            ##################
##############################################################################
##############################################################################

DEGs.analysis.specific.comparison <- function(contrasts.results,
                                              number.comparison, 
                                              threshold.log.fc,
                                              threshold.q.value,
                                              file.name.base)
{
  
  volcano.file.name <- paste(c(file.name.base,"_volcano",".png"),collapse = "")
  file.name.activated.genes <- paste(c(file.name.base,"_activated",".tsv"),collapse = "")
  file.name.repressed.genes <- paste(c(file.name.base,"_repressed",".tsv"),collapse = "")
  barplot.file.name.act <- paste(c(file.name.base,"_activated_barplot",".png"),collapse = "")
  dotplot.file.name.act <- paste(c(file.name.base,"_activated_dotplot",".png"),collapse = "")
  cnetplot.file.name.act <- paste(c(file.name.base,"_activated_cnetplot",".png"),collapse = "")
  barplot.file.name.rep <- paste(c(file.name.base,"_repressed_barplot",".png"),collapse = "")
  dotplot.file.name.rep <- paste(c(file.name.base,"_repressed_dotplot",".png"),collapse = "")
  cnetplot.file.name.rep <- paste(c(file.name.base,"_repressed_cnetplot",".png"),collapse = "")
  GO.file.name.activated.genes <- paste(c(file.name.base,"_activated_GO",".tsv"),collapse = "")
  GO.file.name.repressed.genes <- paste(c(file.name.base,"_repressed_GO",".tsv"),collapse = "")
  
  
  #Sacamos los datos del primer contraste en un nuevo objeto
  comparison.results = topTable(contrasts.results, number = 33602, 
                                coef = number.comparison, sort.by = "logFC")
  
  #sacamos su fold change
  comparison.results.logFC = comparison.results$logFC
  #sacamos su q-value
  comparison.results.qValue = comparison.results$adj.P.Val
  
  #sacamos los ids de los genes y se los añadimos al FC y al q-value
  comparison.results.gene.ids = rownames(x = comparison.results)
  names(x = comparison.results.logFC) = comparison.results.gene.ids
  names(x = comparison.results.qValue) = comparison.results.gene.ids
  #Identificamos los genes activados 
  comparison.activated.genes = comparison.results.gene.ids[comparison.results.logFC > threshold.log.fc & 
                                                             comparison.results.qValue < threshold.q.value]
  length(x = comparison.activated.genes) ## Con logFC > 1 y q-value > 0.05 salen: 777 genes activados
  #Identificamos los genes reprimidos 
  comparison.repressed.genes = comparison.results.gene.ids[comparison.results.logFC < (-threshold.log.fc) & 
                                                             comparison.results.qValue < threshold.q.value]
  length(x = comparison.repressed.genes) ## Con logFC < -1 y q-value > 0.05 salen: 463 gener reprimidos
  ## Volcano plot
  comparison.results.logQval = -log10(x = comparison.results.qValue)
  
  png(filename = volcano.file.name)
  plot(comparison.results.logFC, comparison.results.logQval, pch = 19, col = "grey", cex = 0.8,
       #xlim = c(-7,7), ylim = c(0,4), 
       xlab = "log2(Fold-chage)", ylab = "-log10(q-value)", cex.lab = 1.5)
  points(x = comparison.results.logFC[comparison.activated.genes], 
         y = comparison.results.logQval[comparison.activated.genes], pch = 19, 
         col = "dark green", cex = 0.8)
  points(x = comparison.results.logFC[comparison.repressed.genes], 
         y = comparison.results.logQval[comparison.repressed.genes], pch = 19, 
         col = "brown", cex = 0.8)
  lines(c(1, 1), c(0,40), col = "black", lwd = 2.5, lty = "dashed")
  lines(c(-1, -1), c(0,40), col = "black", lwd = 2.5, lty = "dashed")
  lines(c(-7:7), rep(x = -log10(x = 0.05), 15), col = "black", lwd = 2.5, lty = "dashed")
  dev.off()
  
  write.table(comparison.activated.genes, file = file.name.activated.genes, 
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  write.table(comparison.repressed.genes, file = file.name.repressed.genes, 
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  
  
  ## Enriquecimiento funcional: 
  
  ## Activated:
  comparison.activated.enrichGO = enrichGO(gene = comparison.activated.genes,
                                           OrgDb = org.At.tair.db, ont = "BP",
                                           pAdjustMethod = "BH", pvalueCutoff = 0.05,
                                           readable = TRUE, keyType = "TAIR")
  write.table(comparison.activated.enrichGO, file =  GO.file.name.activated.genes,
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  
  
  #SI PONES READABLE = TRUE, PUEDES HACER UN WRITE.TABLE Y SE GUARDA TAMBIEN PARA 
  #CADA GOENRICHMENT
  png(filename = barplot.file.name.act)
  print(barplot(comparison.activated.enrichGO, showCategory = 25))
  dev.off()
  
  png(filename = dotplot.file.name.act)
  print(dotplot(comparison.activated.enrichGO, showCategory = 20))
  dev.off()
  png(filename = cnetplot.file.name.act)
  print(cnetplot(comparison.activated.enrichGO, showCategory = 20))
  dev.off()
  
  ## Repressed:
  comparison.repressed.enrichGO = enrichGO(gene = comparison.repressed.genes, 
                                           OrgDb = org.At.tair.db, ont = "BP",
                                           pAdjustMethod = "BH", pvalueCutoff = 0.05,
                                           readable = FALSE, keyType = "TAIR")
  
  png(filename = barplot.file.name.rep)
  print(barplot(comparison.repressed.enrichGO, showCategory = 25))
  dev.off()
  png(filename = dotplot.file.name.rep)
  print(dotplot(comparison.repressed.enrichGO, showCategory = 20))
  dev.off()
  png(filename = cnetplot.file.name.rep)
  print(cnetplot(comparison.repressed.enrichGO, showCategory = 20))
  dev.off()
  
  
  ## Enriquecimiento rutas metabólicas (KEGG): 
  comparison.activated.enrichKEGG = enrichKEGG(gene = comparison.activated.genes,
                                               organism = "ath",
                                               pAdjustMethod = "BH",
                                               pvalueCutoff  = 0.05)
  comparison.activated.enrichKEGG.df = as.data.frame(comparison.activated.enrichKEGG)
  
  
  
  comparison.repressed.enrichKEGG = enrichKEGG(gene = comparison.repressed.genes,
                                               organism = "ath",
                                               pAdjustMethod = "BH",
                                               pvalueCutoff = 0.05)
  comparison.repressed.enrichKEGG.df = as.data.frame(comparison.repressed.enrichKEGG)
  
  return(list(activated=comparison.activated.genes,repressed=comparison.repressed.genes))
}


##############################################################################
##############################################################################
################                                            ##################
################               FIN DE LA FUNCIÓN            ##################
################                                            ##################
##############################################################################
##############################################################################


##El paquete **limma** (Linear Models for Microarray Analysis) proporciona las 
##funciones necesarias para determinar los genes expresados de forma 
##diferencial (DEGs). 
library(limma)

## Especificamos el diseño experimental

experimental.design = model.matrix(~ -1+factor(c(1,1,1,2,2,2,3,3,3,4,4,4,
                                                 5,5,6,6,6,7,7,7,8,8,8,
                                                 9,9,9,10,10,10,11,11,11,12,12,12,
                                                 13,13,13,14,14,14,15,15,15,16,16,16,
                                                 17,17,17,18,18,18,19,19,19,20,20,20,
                                                 21,21,21,22,22,22,23,23,23,24,24,24,
                                                 25,25,25,26,26,26,27,27,27,28,28,28)))
colnames(experimental.design) = c("ex1velm", "ex1velt", "ex2velm",  "ex2velt",
                                  "ex3velm", "ex3velt",  "ex4velm", "ex4velt",
                                  "ex1bonm", "ex1bont", "ex2bonm", "ex2bont",
                                  "ex3bonm", "ex3bont", "ex1mojm",  "ex1mojt",
                                  "ex2mojm", "ex2mojt",  "ex3mojm",  "ex3mojt",
                                  "ex1col0m","ex1col0t", "ex2col0m", "ex2col0t",
                                  "ex3col0m", "ex3col0t", "ex4col0m", "ex4col0t")


# A continuación, ajustamos la estimación de los niveles de expresión de cada
## gen a un modelo lineal teniendo en cuenta el diseño experimental. Este paso
## fundamentalmente se calcula la media de las réplicas en cada condición.

linear.fit = lmFit(log.normalized.ecofun.gene.expression, experimental.design)

##Para especificar los constrastes a realizar utilizamos la función
##*makeContrasts* que recibe como entrada los contrastes a realizar separados 
##por comas y especificados con los nombres de las dos condiciones 
##correspondientes separadas por un guión -. También recibe el argumento 
##levels, un vector con el nombre de las condiciones:


#Primera tanda de contrastes: tarde vs. mañana(control) 
#(dentro de cada ecotipo y de cada muestreo)
contrast.matrix <- makeContrasts(ex1velt - ex1velm, ex2velt - ex2velm,
                                 ex3velt - ex3velm, ex4velt - ex4velm, 
                                 ex1bont - ex1bonm, ex2bont - ex2bonm, 
                                 ex3bont - ex3bonm, ex1mojt - ex1mojm, 
                                 ex2mojt - ex2mojm, ex3mojt - ex3mojm,
                                 ex1col0t - ex1col0m, ex2col0t - ex2col0m, 
                                 ex3col0t - ex3col0m, ex4col0t - ex4col0m, 
                                 levels=c("ex1velm", "ex1velt","ex2velm", "ex2velt",
                                          "ex3velm", "ex3velt","ex4velm", "ex4velt",
                                          "ex1bonm", "ex1bont","ex2bonm", "ex2bont",
                                          "ex3bonm", "ex3bont", "ex1mojm", "ex1mojt",
                                          "ex2mojm", "ex2mojt","ex3mojm", "ex3mojt", 
                                          "ex1col0m", "ex1col0t","ex2col0m", "ex2col0t",
                                          "ex3col0m", "ex3col0t","ex4col0m", "ex4col0t"))

##Calculamos el fold-change y los p-valores correspondientes para cada gen en
##cada uno de los constrastes especificados utilizando las funciones *constrasts.fit* 
##y *eBayes*.

contrasts.linear.fit = contrasts.fit(linear.fit, contrasts = contrast.matrix)
contrasts.results = eBayes(contrasts.linear.fit)
nrow(x = log.normalized.ecofun.gene.expression)

#contrastes individuales tarde vs. mañana:
degs.ex1velt_m <- DEGs.analysis.specific.comparison(contrasts.results = contrasts.results,
                                                    number.comparison = 1, 
                                                    threshold.log.fc = 1,
                                                    threshold.q.value = 0.05,
                                                    file.name.base = "results/ex1velt_m")

DEGs.analysis.specific.comparison(contrasts.results = contrasts.results,
                                  number.comparison = 2, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex2velt_m")

DEGs.analysis.specific.comparison(contrasts.results = contrasts.results,
                                  number.comparison = 3, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex3velt_m")

DEGs.analysis.specific.comparison(contrasts.results = contrasts.results,
                                  number.comparison = 4, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex4velt_m")

DEGs.analysis.specific.comparison(contrasts.results = contrasts.results,
                                  number.comparison = 5, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex1bont_m")

DEGs.analysis.specific.comparison(contrasts.results = contrasts.results,
                                  number.comparison = 6, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex2bont_m")

DEGs.analysis.specific.comparison(contrasts.results = contrasts.results,
                                  number.comparison = 7, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex3bont_m")

DEGs.analysis.specific.comparison(contrasts.results = contrasts.results,
                                  number.comparison = 8, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex1mojt_m")

DEGs.analysis.specific.comparison(contrasts.results = contrasts.results,
                                  number.comparison = 9, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex2mojt_m")

DEGs.analysis.specific.comparison(contrasts.results = contrasts.results,
                                  number.comparison = 10, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex3mojt_m")

DEGs.analysis.specific.comparison(contrasts.results = contrasts.results,
                                  number.comparison = 11, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex1col0t_m")

DEGs.analysis.specific.comparison(contrasts.results = contrasts.results,
                                  number.comparison = 12, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex2col0t_m")

DEGs.analysis.specific.comparison(contrasts.results = contrasts.results,
                                  number.comparison = 13, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex3col0t_m")

DEGs.analysis.specific.comparison(contrasts.results = contrasts.results,
                                  number.comparison = 14, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex4col0t_m")



#Segundaa tanda de contrastes: dentro de cada ecotipo entre cada muestreo, 
#el muestreo anterior será el control
contrast.matrix2 <- makeContrasts(ex2velm - ex1velm, ex2velt - ex1velt, ex3velm - ex1velm,
                                  ex3velt - ex1velt, ex3velm - ex2velm, ex3velt - ex2velt,
                                  ex4velm - ex1velm, ex4velt - ex1velt, ex4velm - ex2velm,
                                  ex4velt - ex2velt, ex4velm - ex3velm, ex4velt - ex3velt,
                                  ex2bonm - ex1bonm, ex2bont - ex1bont, ex3bonm - ex1bonm,
                                  ex3bont - ex1bont, ex3bonm - ex2bonm, ex3bont - ex2bont,
                                  ex2mojm - ex1mojm, ex2mojt - ex1mojt, ex3mojm - ex1mojm,
                                  ex3mojt - ex1mojt, ex3mojm - ex2mojm, ex3mojt - ex2mojt,
                                  ex2col0m - ex1col0m, ex2col0t - ex1col0t, ex3col0m - ex1col0m,
                                  ex3col0t - ex1col0t, ex3col0m - ex2col0m, ex3col0t - ex2col0t,
                                  ex4col0m - ex1col0m, ex4col0t - ex1col0t, ex4col0m - ex2col0m,
                                  ex4col0t - ex2col0t, ex4col0m - ex3col0m, ex4col0t - ex3col0t,
                                  levels=c("ex1velm", "ex1velt","ex2velm", "ex2velt",
                                           "ex3velm", "ex3velt","ex4velm", "ex4velt",
                                           "ex1bonm", "ex1bont","ex2bonm", "ex2bont",
                                           "ex3bonm", "ex3bont", "ex1mojm", "ex1mojt",
                                           "ex2mojm", "ex2mojt","ex3mojm", "ex3mojt", 
                                           "ex1col0m", "ex1col0t","ex2col0m", "ex2col0t",
                                           "ex3col0m", "ex3col0t","ex4col0m", "ex4col0t"))

##Calculamos el fold-change y los p-valores correspondientes para cada gen en
##cada uno de los constrastes especificados utilizando las funciones *constrasts.fit* 
##y *eBayes*.

contrasts.linear.fit2 = contrasts.fit(linear.fit, contrasts = contrast.matrix2)
contrasts.results2 = eBayes(contrasts.linear.fit2)
nrow(x = log.normalized.ecofun.gene.expression)

#Contrastes para veleta entre muestreos:
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results2,
                                  number.comparison = 1, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex2-ex1velm")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results2,
                                  number.comparison = 2, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex2-ex1velt")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results2,
                                  number.comparison = 3, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex3-ex1velm")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results2,
                                  number.comparison = 4, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex3-ex1velt")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results2,
                                  number.comparison = 5, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex3-ex2velm")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results2,
                                  number.comparison = 6, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.1,
                                  file.name.base = "results/ex3-ex2velt")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results2,
                                  number.comparison = 7, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex4-ex1velm")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results2,
                                  number.comparison = 8, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex4-ex1velt")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results2,
                                  number.comparison = 9, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex4-ex2velm")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results2,
                                  number.comparison = 10, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex4-ex2velt")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results2,
                                  number.comparison = 11, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex4-ex3velm")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results2,
                                  number.comparison = 12, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex4-ex3velt")

#Contrastes para Bonanza entre muestreos:
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results2,
                                  number.comparison = 13, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex2-ex1bonm")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results2,
                                  number.comparison = 14, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex2-ex1bont")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results2,
                                  number.comparison = 15, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex3-ex1bonm")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results2,
                                  number.comparison = 16, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex3-ex1bont")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results2,
                                  number.comparison = 17, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex3-ex2bonm")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results2,
                                  number.comparison = 18, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex3-ex2bont")

#Contrastes para Montejaque entre muestreos:
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results2,
                                  number.comparison = 19, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex2-ex1mojm")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results2,
                                  number.comparison = 20, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex2-ex1mojt")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results2,
                                  number.comparison = 21, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex3-ex1mojm")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results2,
                                  number.comparison = 22, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex3-ex1mojt")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results2,
                                  number.comparison = 23, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex3-ex2mojm")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results2,
                                  number.comparison = 24, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex3-ex2mojt")

#Contrastes para veleta entre muestreos:
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results2,
                                  number.comparison = 25, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex2-ex1col0m")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results2,
                                  number.comparison = 26, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex2-ex1col0t")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results2,
                                  number.comparison = 27, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex3-ex1col0m")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results2,
                                  number.comparison = 28, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex3-ex1col0t")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results2,
                                  number.comparison = 29, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex3-ex2col0m")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results2,
                                  number.comparison = 30, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex3-ex2col0t")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results2,
                                  number.comparison = 31, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex4-ex1col0m")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results2,
                                  number.comparison = 32, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex4-ex1col0t")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results2,
                                  number.comparison = 33, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex4-ex2col0m")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results2,
                                  number.comparison = 34, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex4-ex2col0t")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results2,
                                  number.comparison = 35, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex4-ex3col0m")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results2,
                                  number.comparison = 36, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex4-ex3col0t")




#Tercera tanda de contrastes: ecotivos vs. Montejaque(control) 
#(para cada muestreo y cada mañana y tarde)
contrast.matrix3 <- makeContrasts(ex1velm - ex1mojm, ex1velt - ex1mojt,
                                  ex2velm - ex2mojm, ex2velt - ex2mojt,
                                  ex3velm - ex3mojm, ex3velt - ex3mojt,
                                  ex1bonm - ex1mojm, ex1bont - ex1mojt,
                                  ex2bonm - ex2mojm, ex2bont - ex2mojt,
                                  ex3bonm - ex3mojm, ex3bont - ex3mojt,
                                  ex1col0m - ex1mojm, ex1col0t - ex1mojt,
                                  ex2col0m - ex2mojm, ex2col0t - ex2mojt,
                                  ex3col0m - ex3mojm, ex3col0t - ex3mojt,
                                  levels=c("ex1velm", "ex1velt","ex2velm", "ex2velt",
                                           "ex3velm", "ex3velt","ex4velm", "ex4velt",
                                           "ex1bonm", "ex1bont","ex2bonm", "ex2bont",
                                           "ex3bonm", "ex3bont", "ex1mojm", "ex1mojt",
                                           "ex2mojm", "ex2mojt","ex3mojm", "ex3mojt", 
                                           "ex1col0m", "ex1col0t","ex2col0m", "ex2col0t",
                                           "ex3col0m", "ex3col0t","ex4col0m", "ex4col0t"))

##Calculamos el fold-change y los p-valores correspondientes para cada gen en
##cada uno de los constrastes especificados utilizando las funciones *constrasts.fit* 
##y *eBayes*.

contrasts.linear.fit3 = contrasts.fit(linear.fit, contrasts = contrast.matrix3)
contrasts.results3 = eBayes(contrasts.linear.fit3)
nrow(x = log.normalized.ecofun.gene.expression)

#Contrastes para veleta vs montejaque:
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results3,
                                  number.comparison = 1, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex1vel-mojm")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results3,
                                  number.comparison = 2, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex1vel-mojt")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results3,
                                  number.comparison = 3, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex2vel-mojm")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results3,
                                  number.comparison = 4, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex2vel-mojt")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results3,
                                  number.comparison = 5, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex3vel-mojm")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results3,
                                  number.comparison = 6, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex3vel-mojt")

#Contrastes para bonanza vs montejaque:
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results3,
                                  number.comparison = 7, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex1bon-mojm")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results3,
                                  number.comparison = 8, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex1bon-mojt")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results3,
                                  number.comparison = 9, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex2bon-mojm")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results3,
                                  number.comparison = 10, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex2bon-mojt")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results3,
                                  number.comparison = 11, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex3bon-mojm")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results3,
                                  number.comparison = 12, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex3bon-mojt")

#Contrastes para columbia vs montejaque:
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results3,
                                  number.comparison = 13, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex1col0-mojm")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results3,
                                  number.comparison = 14, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex1col0-mojt")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results3,
                                  number.comparison = 15, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex2col0-mojm")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results3,
                                  number.comparison = 16, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex2col0-mojt")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results3,
                                  number.comparison = 17, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex3col0-mojm")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results3,
                                  number.comparison = 18, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex3col0-mojt")


#Cuarta tanda de contrastes: ecotivos vs. Columbia(control) 
#(para cada muestreo y cada mañana y tarde)
contrast.matrix4 <- makeContrasts(ex1velm - ex1col0m, ex1velt - ex1col0t,
                                  ex2velm - ex2col0m, ex2velt - ex2col0t,
                                  ex3velm - ex3col0m, ex3velt - ex3col0t,
                                  ex4velt - ex4col0t,
                                  ex1bonm - ex1col0m, ex1bont - ex1col0t,
                                  ex2bonm - ex2col0m, ex2bont - ex2col0t,
                                  ex3bonm - ex3col0m, ex3bont - ex3col0t,
                                  ex1mojm - ex1col0m, ex1mojt - ex1col0t,
                                  ex2mojm - ex2col0m, ex2mojt - ex2col0t,
                                  ex3mojm - ex3col0m, ex3mojt - ex3col0t,
                                  levels=c("ex1velm", "ex1velt","ex2velm", "ex2velt",
                                           "ex3velm", "ex3velt","ex4velm", "ex4velt",
                                           "ex1bonm", "ex1bont","ex2bonm", "ex2bont",
                                           "ex3bonm", "ex3bont", "ex1mojm", "ex1mojt",
                                           "ex2mojm", "ex2mojt","ex3mojm", "ex3mojt", 
                                           "ex1col0m", "ex1col0t","ex2col0m", "ex2col0t",
                                           "ex3col0m", "ex3col0t","ex4col0m", "ex4col0t"))

##Calculamos el fold-change y los p-valores correspondientes para cada gen en
##cada uno de los constrastes especificados utilizando las funciones *constrasts.fit* 
##y *eBayes*.

contrasts.linear.fit4 = contrasts.fit(linear.fit, contrasts = contrast.matrix4)
contrasts.results4 = eBayes(contrasts.linear.fit4)
nrow(x = log.normalized.ecofun.gene.expression)

#Contrastes para veleta vs columbia:
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results4,
                                  number.comparison = 1, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex1vel-col0m")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results4,
                                  number.comparison = 2, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex1vel-col0t")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results4,
                                  number.comparison = 3, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex2vel-col0m")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results4,
                                  number.comparison = 4, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex2vel-col0t")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results4,
                                  number.comparison = 5, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex3vel-col0m")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results4,
                                  number.comparison = 6, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex3vel-col0t")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results4,
                                  number.comparison = 7, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex4vel-col0t")


#Contrastes para bonanza vs columbia:
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results4,
                                  number.comparison = 8, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex1bon-col0m")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results4,
                                  number.comparison = 9, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex1bon-col0t")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results4,
                                  number.comparison = 10, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex2bon-col0m")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results4,
                                  number.comparison = 11, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex2bon-col0t")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results4,
                                  number.comparison = 12, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex3bon-col0m")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results4,
                                  number.comparison = 13, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex3bon-col0t")

#Contrastes para montejaque vs columbia:
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results4,
                                  number.comparison = 14, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex1moj-col0m")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results4,
                                  number.comparison = 15, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex1moj-col0t")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results4,
                                  number.comparison = 16, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex2moj-col0m")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results4,
                                  number.comparison = 17, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex2moj-col0t")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results4,
                                  number.comparison = 18, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex3moj-col0m")
DEGs.analysis.specific.comparison(contrasts.results = contrasts.results4,
                                  number.comparison = 19, 
                                  threshold.log.fc = 1,
                                  threshold.q.value = 0.05,
                                  file.name.base = "results/ex3moj-col0t")



##PRIMERA TANDA DE CONTRASTES (TARDE VS MAÑANA(CONTROL))
#ENTRE ECOTIPOS:
## Activados:
library(VennDiagram)
dev.off
grid.newpage()
activated.venn = venn.diagram(x = list(col0 = read.table("results/ex1col0t_m_activated.tsv", 
                                                         as.is = TRUE, header = FALSE)[[1]], 
                                       moj = read.table("results/ex1mojt_m_activated.tsv", 
                                                        as.is = TRUE, header = FALSE)[[1]],
                                       vel = read.table("results/ex1velt_m_activated.tsv", 
                                                        as.is = TRUE, header = FALSE)[[1]],
                                       bon = read.table("results/ex1bont_m_activated.tsv", 
                                                        as.is = TRUE, header = FALSE)[[1]]),
                              category.names = c("Col0", "Moj", "Vel","Bon"), output = TRUE,
                              filename = NULL, lwd = 2, lty = 'blank', 
                              fill = c("light green", "light blue", "red", "yellow"))
grid.draw(x = activated.venn)


## Inhibidos:
grid.newpage()
repressed.venn = venn.diagram(x = list(col0 = read.table("results/ex1col0t_m_repressed.tsv",
                                                         as.is = TRUE, header = FALSE)[[1]], 
                                       moj = read.table("results/ex1mojt_m_repressed.tsv",
                                                        as.is = TRUE, header = FALSE)[[1]],
                                       vel = read.table("results/ex1velt_m_repressed.tsv",
                                                        as.is = TRUE, header = FALSE)[[1]],
                                       bon = read.table("results/ex1bont_m_repressed.tsv",
                                                        as.is = TRUE, header = FALSE)[[1]]), 
                              category.names = c("Col0", "Moj", "Vel","Bon"), output = TRUE,
                              filename = NULL, lwd = 2, lty = 'blank', 
                              fill = c("light green", "light blue", "red", "yellow")) 

grid.draw(x = repressed.venn)

#common genes

col0 = read.table("results/ex1col0t_m_activated.tsv", 
                  as.is = TRUE, header = FALSE)[[1]]
moj = read.table("results/ex1mojt_m_activated.tsv", 
                 as.is = TRUE, header = FALSE)[[1]]
vel = read.table("results/ex1velt_m_activated.tsv", 
                 as.is = TRUE, header = FALSE)[[1]]
bon = read.table("results/ex1bont_m_activated.tsv", 
                 as.is = TRUE, header = FALSE)[[1]]

inter = intersect(intersect(intersect(col0, moj), vel),bon)
length(inter)
vel.espc = setdiff(vel, unique(c(col0,moj,bon)))
length(vel.espc)
#hacer enrich.go de vel.espec para ver qué funciones tienen los genes 
#especificos de veleta y cuales son estos genes

write(x = vel.espc, file = "results/genes_ex1t_m_vel.tsv")



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # FUNCIÓN DIAGRAMAS DE VENN, GENES COMUNES Y ENRIQUECIMIENTO GO # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
venn.cg.enrgo = function(grupo, group.names, nombre)
{
  #nombres diagramas de venn
  venn.file.name.act = paste(c("results/", grupo,"v_",nombre,"_act",".png"),
                             collapse = "")
  venn.file.name.rep = paste(c("results/", grupo,"v_",nombre,"_rep",".png"),
                             collapse = "")
  #nombres genes comunes y específicos
  cg.file.name.act = paste(c("results/", grupo,"cg_",nombre,"_act",".tsv"),
                           collapse = "")
  eg1.file.name.act = paste(c("results/", grupo,"eg_",group.names[1],"_act",".tsv"),
                            collapse = "")
  eg2.file.name.act = paste(c("results/", grupo,"eg_",group.names[2],"_act",".tsv"),
                            collapse = "")
  eg3.file.name.act = paste(c("results/", grupo,"eg_",group.names[3],"_act",".tsv"),
                            collapse = "")
  eg4.file.name.act =  paste(c("results/",grupo,"eg_",group.names[4],"_act",".tsv"),
                             collapse = "")
  cg.file.name.rep = paste(c("results/",grupo,"cg_",nombre,"_rep",".tsv"),
                           collapse = "")
  eg1.file.name.rep = paste(c("results/",grupo,"eg_",group.names[1],"_rep",".tsv"),
                            collapse = "")
  eg2.file.name.rep = paste(c("results/",grupo,"eg_",group.names[2],"_rep",".tsv"),
                            collapse = "")
  eg3.file.name.rep = paste(c("results/",grupo,"eg_",group.names[3],"_rep",".tsv"),
                            collapse = "")
  eg4.file.name.rep = paste(c("results/",grupo,"eg_",group.names[4],"_rep",".tsv"),
                            collapse = "")
  #nombres enriquecimiento funcional
  cg.file.name.actGO = paste(c("results/",grupo,"cg_",nombre,"_actGO",".tsv"),
                             collapse = "")
  eg1.file.name.actGO = paste(c("results/",grupo,"eg_",group.names[1],"_actGO",".tsv"),
                              collapse = "")
  eg2.file.name.actGO = paste(c("results/",grupo,"eg_",group.names[2],"_actGO",".tsv"),
                              collapse = "")
  eg3.file.name.actGO = paste(c("results/",grupo,"eg_",group.names[3],"_actGO",".tsv"),
                              collapse = "")
  eg4.file.name.actGO = paste(c("results/",grupo,"eg_",group.names[4],"_actGO",".tsv"),
                              collapse = "")
  cg.file.name.repGO = paste(c("results/",grupo,"cg_",nombre,"_repGO",".tsv"),
                             collapse = "")
  eg1.file.name.repGO = paste(c("results/",grupo,"eg_",group.names[1],"_repGO",".tsv"),
                              collapse = "")
  eg2.file.name.repGO = paste(c("results/",grupo,"eg_",group.names[2],"_repGO",".tsv"),
                              collapse = "")
  eg3.file.name.repGO = paste(c("results/",grupo,"eg_",group.names[3],"_repGO",".tsv"),
                              collapse = "")
  eg4.file.name.repGO = paste(c("results/",grupo,"eg_",group.names[4],"_repGO",".tsv"),
                              collapse = "")
  #nombres dde los barplots
  cg.file.name.actBP = paste(c("results/",grupo,"cg_",nombre,"_act_bp",".png"),
                             collapse = "")
  eg1.file.name.actBP = paste(c("results/",grupo,"eg_",group.names[1],"_act_bp",".png"),
                              collapse = "")
  eg2.file.name.actBP = paste(c("results/",grupo,"eg_",group.names[2],"_act_bp",".png"),
                              collapse = "")
  eg3.file.name.actBP = paste(c("results/",grupo,"eg_",group.names[3],"_act_bp",".png"),
                              collapse = "")
  eg4.file.name.actBP = paste(c("results/",grupo,"eg_",group.names[4],"_act_bp",".png"),
                              collapse = "")
  cg.file.name.repBP = paste(c("results/",grupo,"cg_",nombre,"_rep_bp",".png"),
                             collapse = "")
  eg1.file.name.repBP = paste(c("results/",grupo,"eg_",group.names[1],"_rep_bp",".png"),
                              collapse = "")
  eg2.file.name.repBP = paste(c("results/",grupo,"eg_",group.names[2],"_rep_bp",".png"),
                              collapse = "")
  eg3.file.name.repBP = paste(c("results/",grupo,"eg_",group.names[3],"_rep_bp",".png"),
                              collapse = "")
  eg4.file.name.repBP = paste(c("results/",grupo,"eg_",group.names[4],"_rep_bp",".png"),
                              collapse = "")
  
  if (length(group.names) == 4) {
    #aquí meto la parte de la función para cuatro muestras
    fst_act = read.table((paste(c("results/",group.names[1],"_activated",".tsv"),collapse = "")), 
                         as.is = TRUE, header = FALSE)[[1]]
    scn_act = read.table((paste(c("results/",group.names[2],"_activated",".tsv"),collapse = "")), 
                         as.is = TRUE, header = FALSE)[[1]]
    thr_act = read.table((paste(c("results/",group.names[3],"_activated",".tsv"),collapse = "")), 
                         as.is = TRUE, header = FALSE)[[1]]
    for_act = read.table((paste(c("results/",group.names[4],"_activated",".tsv"),collapse = "")), 
                         as.is = TRUE, header = FALSE)[[1]]
    fst_rep = read.table((paste(c("results/",group.names[1],"_repressed",".tsv"),collapse = "")), 
                         as.is = TRUE, header = FALSE)[[1]]
    scn_rep = read.table((paste(c("results/",group.names[2],"_repressed",".tsv"),collapse = "")), 
                         as.is = TRUE, header = FALSE)[[1]]
    thr_rep = read.table((paste(c("results/",group.names[3],"_repressed",".tsv"),collapse = "")), 
                         as.is = TRUE, header = FALSE)[[1]]
    for_rep = read.table((paste(c("results/",group.names[4],"_repressed",".tsv"),collapse = "")), 
                         as.is = TRUE, header = FALSE)[[1]]
    #VENN ACTIVADOS
    png(filename = venn.file.name.act, width = 500, height = 500)
    grid.newpage()
    activated.venn = venn.diagram(x = list(fst_act, scn_act, thr_act, for_act),
                                  category.names = group.names, output = TRUE,
                                  filename = NULL, lwd = 2, lty = 'blank', 
                                  fill = c("light green", "light blue", "red", "yellow"))
    grid.draw(x = activated.venn)
    dev.off()
    ## Inhibidos:
    png(filename = venn.file.name.rep, width = 500, height = 500)
    grid.newpage()
    repressed.venn = venn.diagram(x = list(fst_rep, scn_rep, thr_rep, for_rep), 
                                  category.names = group.names, output = TRUE,
                                  filename = NULL, lwd = 2, lty = 'blank', 
                                  fill = c("light green", "light blue", "red", "yellow")) 
    grid.draw(x = repressed.venn)
    dev.off()
    
    ##genes comunes y específicos:
    genes_comunes = intersect(intersect(intersect(fst_act, scn_act), thr_act),for_act)
    write(x = genes_comunes, file = cg.file.name.act)
    
    fst.espc = setdiff(fst_act, unique(c(scn_act,thr_act,for_act)))
    write(x = fst.espc, file = eg1.file.name.act)
    
    scn.espc = setdiff(scn_act, unique(c(fst_act,thr_act,for_act)))
    write(x = fst.espc, file = eg2.file.name.act)
    
    thr.espc = setdiff(thr_act, unique(c(scn_act,fst_act,for_act)))
    write(x = fst.espc, file = eg3.file.name.act)
    
    for.espc = setdiff(for_act, unique(c(scn_act,thr_act,fst_act)))
    write(x = fst.espc, file = eg4.file.name.act)
    
    ## Enriquecimiento funcional: 
    ## Common Activated:
    comparison.activated.enrichGO = enrichGO(gene = genes_comunes,
                                             OrgDb = org.At.tair.db, ont = "BP",
                                             pAdjustMethod = "BH", pvalueCutoff = 0.05,
                                             readable = TRUE, keyType = "TAIR")
    write.table(comparison.activated.enrichGO,file = cg.file.name.actGO ,
                quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
    png(filename = cg.file.name.actBP)
    print(barplot(comparison.activated.enrichGO, showCategory = 25))
    dev.off()
    ## specific Activated:
    #primero:
    comparison.activated.enrichGO2 = enrichGO(gene = fst.espc,
                                              OrgDb = org.At.tair.db, ont = "BP",
                                              pAdjustMethod = "BH", pvalueCutoff = 0.05,
                                              readable = TRUE, keyType = "TAIR")
    write.table(comparison.activated.enrichGO2,file = eg1.file.name.act,
                quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
    png(filename = eg1.file.name.actBP)
    print(barplot(comparison.activated.enrichGO2, showCategory = 25))
    dev.off()
    #segundo:
    comparison.activated.enrichGO3 = enrichGO(gene = scn.espc,
                                              OrgDb = org.At.tair.db, ont = "BP",
                                              pAdjustMethod = "BH", pvalueCutoff = 0.05,
                                              readable = TRUE, keyType = "TAIR")
    write.table(comparison.activated.enrichGO3,file =  eg2.file.name.actGO,
                quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
    png(filename = eg2.file.name.actBP)
    print(barplot(comparison.activated.enrichGO3, showCategory = 25))
    dev.off()
    #tercero:
    comparison.activated.enrichGO4 = enrichGO(gene = thr.espc,
                                              OrgDb = org.At.tair.db, ont = "BP",
                                              pAdjustMethod = "BH", pvalueCutoff = 0.05,
                                              readable = TRUE, keyType = "TAIR")
    write.table(comparison.activated.enrichGO4,file =  eg3.file.name.actGO,
                quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
    png(filename = eg3.file.name.actBP)
    print(barplot(comparison.activated.enrichGO4, showCategory = 25))
    dev.off()
    #cuarto:
    comparison.activated.enrichGO5 = enrichGO(gene = for.espc,
                                              OrgDb = org.At.tair.db, ont = "BP",
                                              pAdjustMethod = "BH", pvalueCutoff = 0.05,
                                              readable = TRUE, keyType = "TAIR")
    write.table(comparison.activated.enrichGO5,file =  eg4.file.name.actGO,
                quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
    png(filename = eg4.file.name.actBP)
    print(barplot(comparison.activated.enrichGO5, showCategory = 25))
    dev.off() 
    
    #REPRIMIDOS
    ##genes comunes y específicos:
    genes_comunes_r = intersect(intersect(intersect(fst_rep, scn_rep), thr_rep),for_rep)
    write(x = genes_comunes_r, file = cg.file.name.rep)
    
    fst.espc_r = setdiff(fst_rep, unique(c(scn_rep,thr_rep,for_rep)))
    write(x = fst.espc_r, file = eg1.file.name.rep)
    
    scn.espc_r = setdiff(scn_rep, unique(c(fst_rep,thr_rep,for_rep)))
    write(x = fst.espc_r, file = eg2.file.name.rep)
    
    thr.espc_r = setdiff(thr_rep, unique(c(scn_rep,fst_rep,for_rep)))
    write(x = fst.espc_r, file = eg3.file.name.rep)
    
    for.espc_r = setdiff(for_rep, unique(c(scn_rep,thr_rep,fst_rep)))
    write(x = fst.espc_r, file = eg4.file.name.rep)
    ##Common repressed:
    comparison.repressed.enrichGO = enrichGO(gene = genes_comunes_r,
                                             OrgDb = org.At.tair.db, ont = "BP",
                                             pAdjustMethod = "BH", pvalueCutoff = 0.05,
                                             readable = TRUE, keyType = "TAIR")
    write.table(comparison.repressed.enrichGO,file =  cg.file.name.repGO,
                quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
    png(filename = cg.file.name.repBP)
    print(barplot(comparison.repressed.enrichGO, showCategory = 25))
    dev.off()
    ## specific repressed:
    #primero:
    comparison.repressed.enrichGO2 = enrichGO(gene = fst.espc_r,
                                              OrgDb = org.At.tair.db, ont = "BP",
                                              pAdjustMethod = "BH", pvalueCutoff = 0.05,
                                              readable = TRUE, keyType = "TAIR")
    write.table(comparison.repressed.enrichGO2,file =  eg1.file.name.repGO,
                quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
    png(filename = eg1.file.name.repBP)
    print(barplot(comparison.repressed.enrichGO2, showCategory = 25))
    dev.off()
    #segundo:
    comparison.repressed.enrichGO3 = enrichGO(gene = scn.espc_r,
                                              OrgDb = org.At.tair.db, ont = "BP",
                                              pAdjustMethod = "BH", pvalueCutoff = 0.05,
                                              readable = TRUE, keyType = "TAIR")
    write.table(comparison.repressed.enrichGO3,file =  eg2.file.name.repGO,
                quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
    png(filename = eg2.file.name.repBP)
    print(barplot(comparison.repressed.enrichGO3, showCategory = 25))
    dev.off()
    #tercero:
    comparison.repressed.enrichGO4 = enrichGO(gene = thr.espc_r,
                                              OrgDb = org.At.tair.db, ont = "BP",
                                              pAdjustMethod = "BH", pvalueCutoff = 0.05,
                                              readable = TRUE, keyType = "TAIR")
    write.table(comparison.repressed.enrichGO4, file =  eg3.file.name.repGO,
                quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
    png(filename = eg3.file.name.repBP)
    print(barplot(comparison.repressed.enrichGO4, showCategory = 25))
    dev.off()
    #cuarto:
    comparison.repressed.enrichGO5 = enrichGO(gene = for.espc_r,
                                              OrgDb = org.At.tair.db, ont = "BP",
                                              pAdjustMethod = "BH", pvalueCutoff = 0.05,
                                              readable = TRUE, keyType = "TAIR")
    write.table(comparison.repressed.enrichGO5, file =  eg4.file.name.repGO,
                quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
    png(filename = eg4.file.name.repBP)
    print(barplot(comparison.repressed.enrichGO5, showCategory = 25))
    dev.off() 
    
  } else {
    if (length(group.names) == 3) {
      #aquí meto la parte de la función para tres muestras
      fst_act = read.table((paste(c("results/",group.names[1],"_activated",".tsv"),
                                  collapse = "")), 
                           as.is = TRUE, header = FALSE)[[1]]
      scn_act = read.table((paste(c("results/",group.names[2],"_activated",".tsv"),
                                  collapse = "")), 
                           as.is = TRUE, header = FALSE)[[1]]
      thr_act = read.table((paste(c("results/",group.names[3],"_activated",".tsv"),
                                  collapse = "")), 
                           as.is = TRUE, header = FALSE)[[1]]
      
      fst_rep = read.table((paste(c("results/",group.names[1],"_repressed",".tsv"),
                                  collapse = "")), 
                           as.is = TRUE, header = FALSE)[[1]]
      scn_rep = read.table((paste(c("results/",group.names[2],"_repressed",".tsv"),
                                  collapse = "")), 
                           as.is = TRUE, header = FALSE)[[1]]
      thr_rep = read.table((paste(c("results/",group.names[3],"_repressed",".tsv"),
                                  collapse = "")), 
                           as.is = TRUE, header = FALSE)[[1]]
      
      #VENN ACTIVADOS
      png(filename = venn.file.name.act, width = 500, height = 500)
      grid.newpage()
      activated.venn = venn.diagram(x = list(fst_act, scn_act, thr_act),
                                    category.names = group.names, output = TRUE,
                                    filename = NULL, lwd = 2, lty = 'blank', 
                                    fill = c("light green", "light blue", "red"))
      grid.draw(x = activated.venn)
      dev.off()
      ## Inhibidos:
      png(filename = venn.file.name.rep, width = 500, height = 500)
      grid.newpage()
      repressed.venn = venn.diagram(x = list(fst_rep, scn_rep, thr_rep), 
                                    category.names = group.names, output = TRUE,
                                    filename = NULL, lwd = 2, lty = 'blank', 
                                    fill = c("light green", "light blue", "red")) 
      
      grid.draw(x = repressed.venn)
      dev.off()
      
      ##genes comunes y específicos:
      genes_comunes = intersect(intersect(fst_act, scn_act), thr_act)
      write(x = genes_comunes, file = cg.file.name.act)
      
      fst.espc = setdiff(fst_act, unique(c(scn_act,thr_act)))
      write(x = fst.espc, file = eg1.file.name.act)
      
      scn.espc = setdiff(scn_act, unique(c(fst_act,thr_act)))
      write(x = fst.espc, file = eg2.file.name.act)
      
      thr.espc = setdiff(thr_act, unique(c(scn_act,fst_act)))
      write(x = fst.espc, file = eg3.file.name.act)
      
      
      ## Enriquecimiento funcional: 
      ## Common Activated:
      comparison.activated.enrichGO = enrichGO(gene = genes_comunes,
                                               OrgDb = org.At.tair.db, ont = "BP",
                                               pAdjustMethod = "BH", pvalueCutoff = 0.05,
                                               readable = TRUE, keyType = "TAIR")
      write.table(comparison.activated.enrichGO,file =  cg.file.name.actGO,
                  quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
      png(filename = cg.file.name.actBP)
      print(barplot(comparison.activated.enrichGO, showCategory = 25))
      dev.off()
      
      ## specific Activated:
      #primero:
      comparison.activated.enrichGO2 = enrichGO(gene = fst.espc,
                                                OrgDb = org.At.tair.db, ont = "BP",
                                                pAdjustMethod = "BH", pvalueCutoff = 0.05,
                                                readable = TRUE, keyType = "TAIR")
      write.table(comparison.activated.enrichGO2,file =  eg1.file.name.actGO,
                  quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
      png(filename = eg1.file.name.actBP)
      print(barplot(comparison.activated.enrichGO2, showCategory = 25))
      dev.off()
      
      #segundo:
      comparison.activated.enrichGO3 = enrichGO(gene = scn.espc,
                                                OrgDb = org.At.tair.db, ont = "BP",
                                                pAdjustMethod = "BH", pvalueCutoff = 0.05,
                                                readable = TRUE, keyType = "TAIR")
      write.table(comparison.activated.enrichGO3,file =  eg2.file.name.actGO,
                  quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
      png(filename = eg2.file.name.actBP)
      print(barplot(comparison.activated.enrichGO3, showCategory = 25))
      dev.off()
      #tercero:
      comparison.activated.enrichGO4 = enrichGO(gene = thr.espc,
                                                OrgDb = org.At.tair.db, ont = "BP",
                                                pAdjustMethod = "BH", pvalueCutoff = 0.05,
                                                readable = TRUE, keyType = "TAIR")
      write.table(comparison.activated.enrichGO4,file =  eg3.file.name.actGO,
                  quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
      png(filename = eg3.file.name.actBP)
      print(barplot(comparison.activated.enrichGO4, showCategory = 25))
      dev.off()
      
      
      ##genes comunes y específicos:
      genes_comunes_r = intersect(intersect(fst_rep, scn_rep), thr_rep)
      write(x = genes_comunes_r, file = cg.file.name.rep)
      
      fst.espc_r = setdiff(fst_rep, unique(c(scn_rep,thr_rep)))
      write(x = fst.espc_r, file = eg1.file.name.rep)
      
      scn.espc_r = setdiff(scn_rep, unique(c(fst_rep,thr_rep)))
      write(x = fst.espc_r, file = eg2.file.name.rep)
      
      thr.espc_r = setdiff(thr_rep, unique(c(scn_rep,fst_rep)))
      write(x = fst.espc_r, file = eg3.file.name.rep)
      
      
      ##Common repressed:
      comparison.repressed.enrichGO = enrichGO(gene = genes_comunes_r,
                                               OrgDb = org.At.tair.db, ont = "BP",
                                               pAdjustMethod = "BH", pvalueCutoff = 0.05,
                                               readable = TRUE, keyType = "TAIR")
      write.table(comparison.repressed.enrichGO,file =  cg.file.name.repGO,
                  quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
      png(filename = cg.file.name.repBP)
      print(barplot(comparison.repressed.enrichGO, showCategory = 25))
      dev.off()
      
      ## specific repressed:
      #primero:
      comparison.repressed.enrichGO2 = enrichGO(gene = fst.espc_r,
                                                OrgDb = org.At.tair.db, ont = "BP",
                                                pAdjustMethod = "BH", pvalueCutoff = 0.05,
                                                readable = TRUE, keyType = "TAIR")
      write.table(comparison.repressed.enrichGO2, file =  eg1.file.name.repGO,
                  quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
      png(filename = eg1.file.name.repBP)
      print(barplot(comparison.repressed.enrichGO2, showCategory = 25))
      dev.off()
      #segundo:
      comparison.repressed.enrichGO3 = enrichGO(gene = scn.espc_r,
                                                OrgDb = org.At.tair.db, ont = "BP",
                                                pAdjustMethod = "BH", pvalueCutoff = 0.05,
                                                readable = TRUE, keyType = "TAIR")
      write.table(comparison.repressed.enrichGO3,file =  eg2.file.name.repGO,
                  quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
      png(filename = eg2.file.name.repBP)
      print(barplot(comparison.repressed.enrichGO3, showCategory = 25))
      dev.off()
      #tercero:
      comparison.repressed.enrichGO4 = enrichGO(gene = thr.espc_r,
                                                OrgDb = org.At.tair.db, ont = "BP",
                                                pAdjustMethod = "BH", pvalueCutoff = 0.05,
                                                readable = TRUE, keyType = "TAIR")
      write.table(comparison.repressed.enrichGO4,file =  eg3.file.name.repGO,
                  quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
      png(filename = eg3.file.name.repBP)
      print(barplot(comparison.repressed.enrichGO4, showCategory = 25))
      dev.off()
      
      
    } else {
      fst_act = read.table((paste(c("results/",group.names[1],"_activated",".tsv"),
                                  collapse = "")), 
                           as.is = TRUE, header = FALSE)[[1]]
      scn_act = read.table((paste(c("results/",group.names[2],"_activated",".tsv"),
                                  collapse = "")), 
                           as.is = TRUE, header = FALSE)[[1]]
      
      fst_rep = read.table((paste(c("results/",group.names[1],"_repressed",".tsv"),
                                  collapse = "")), 
                           as.is = TRUE, header = FALSE)[[1]]
      scn_rep = read.table((paste(c("results/",group.names[2],"_repressed",".tsv"),
                                  collapse = "")), 
                           as.is = TRUE, header = FALSE)[[1]]
      
      #VENN ACTIVADOS
      png(filename = venn.file.name.act, width = 500, height = 500)
      grid.newpage()
      activated.venn = venn.diagram(x = list(fst_act, scn_act),
                                    category.names = group.names, output = TRUE,
                                    filename = NULL, lwd = 2, lty = 'blank', 
                                    fill = c("light green", "light blue"))
      grid.draw(x = activated.venn, width = 500, height = 500)
      dev.off()
      ## Inhibidos:
      png(filename = venn.file.name.rep)
      grid.newpage()
      repressed.venn = venn.diagram(x = list(fst_rep, scn_rep), 
                                    category.names = group.names, output = TRUE,
                                    filename = NULL, lwd = 2, lty = 'blank', 
                                    fill = c("light green", "light blue")) 
      
      grid.draw(x = repressed.venn)
      dev.off()
      
      #ACTIVADOS:
      ##genes comunes y específicos:
      genes_comunes = intersect(fst_act, scn_act)
      write(x = genes_comunes, file = cg.file.name.act)
      
      fst.espc = setdiff(fst_act, scn_act)
      write(x = fst.espc, file = eg1.file.name.act)
      
      scn.espc = setdiff(scn_act, fst_act)
      write(x = fst.espc, file = eg2.file.name.act)
      
      ## Enriquecimiento funcional: 
      ## Common Activated:
      comparison.activated.enrichGO = enrichGO(gene = genes_comunes,
                                               OrgDb = org.At.tair.db, ont = "BP",
                                               pAdjustMethod = "BH", pvalueCutoff = 0.05,
                                               readable = TRUE, keyType = "TAIR")
      write.table(comparison.activated.enrichGO, file =  cg.file.name.actGO,
                  quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
      png(filename = cg.file.name.actBP)
      print(barplot(comparison.activated.enrichGO, showCategory = 25))
      dev.off()
      
      ## specific Activated:
      #primero:
      comparison.activated.enrichGO2 = enrichGO(gene = fst.espc,
                                                OrgDb = org.At.tair.db, ont = "BP",
                                                pAdjustMethod = "BH", pvalueCutoff = 0.05,
                                                readable = TRUE, keyType = "TAIR")
      write.table(comparison.activated.enrichGO2,file =  eg1.file.name.actGO,
                  quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
      png(filename = eg1.file.name.actBP)
      print(barplot(comparison.activated.enrichGO2, showCategory = 25))
      dev.off()
      #segundo:
      comparison.activated.enrichGO3 = enrichGO(gene = scn.espc,
                                                OrgDb = org.At.tair.db, ont = "BP",
                                                pAdjustMethod = "BH", pvalueCutoff = 0.05,
                                                readable = TRUE, keyType = "TAIR")
      write.table(comparison.activated.enrichGO3,file =  eg2.file.name.actGO,
                  quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
      png(filename = eg2.file.name.actBP)
      print(barplot(comparison.activated.enrichGO3, showCategory = 25))
      dev.off()
      
      #REPRIMIDOS:
      ##genes comunes y específicos:
      genes_comunes_r = intersect(fst_rep, scn_rep)
      write(x = genes_comunes_r, file = cg.file.name.rep)
      
      fst.espc_r = setdiff(fst_rep, scn_rep)
      write(x = fst.espc_r, file = eg1.file.name.rep)
      
      scn.espc_r = setdiff(scn_rep, fst_rep)
      write(x = fst.espc_r, file = eg2.file.name.rep)
      
      
      ##Common repressed:
      comparison.repressed.enrichGO = enrichGO(gene = genes_comunes_r,
                                               OrgDb = org.At.tair.db, ont = "BP",
                                               pAdjustMethod = "BH", pvalueCutoff = 0.05,
                                               readable = TRUE, keyType = "TAIR")
      write.table(comparison.repressed.enrichGO,file =  cg.file.name.repGO,
                  quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
      png(filename = cg.file.name.repBP)
      print(barplot(comparison.repressed.enrichGO, showCategory = 25))
      dev.off()
      
      ## specific repressed:
      #primero:
      comparison.repressed.enrichGO2 = enrichGO(gene = fst.espc_r,
                                                OrgDb = org.At.tair.db, ont = "BP",
                                                pAdjustMethod = "BH", pvalueCutoff = 0.05,
                                                readable = TRUE, keyType = "TAIR")
      write.table(comparison.repressed.enrichGO2,file =  eg1.file.name.repGO,
                  quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
      png(filename = eg1.file.name.repBP)
      print(barplot(comparison.repressed.enrichGO2, showCategory = 25))
      dev.off()
      
      #segundo:
      comparison.repressed.enrichGO3 = enrichGO(gene = scn.espc_r,
                                                OrgDb = org.At.tair.db, ont = "BP",
                                                pAdjustMethod = "BH", pvalueCutoff = 0.05,
                                                readable = TRUE, keyType = "TAIR")
      write.table(comparison.repressed.enrichGO3,file =  eg2.file.name.repGO,
                  quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
      png(filename = eg2.file.name.repBP)
      print(barplot(comparison.repressed.enrichGO3, showCategory = 25))
      dev.off()
    }
  }
  gc()
}

##############################################################################
##############################################################################
################                                            ##################
################               FIN DE LA FUNCIÓN            ##################
################                                            ##################
##############################################################################
##############################################################################


#Corremos la función para los distintos grupos de genes que queremos hacerles 
#el diagrama de venn y el enriquecimiento:


#TARDE VS MAÑANA:
#GRUPO 1 (g1/): tardes vs mañana entre ECOTIPOS y por cada extracción
tarde_vs_mañana_ex1 = list("ex1col0t_m", "ex1mojt_m", "ex1velt_m", "ex1bont_m")
venn.cg.enrgo("g1/",tarde_vs_mañana_ex1, "t_m_ex1")

tarde_vs_mañana_ex2 = list("ex2col0t_m", "ex2mojt_m", "ex2velt_m", "ex2bont_m")
venn.cg.enrgo("g1/",tarde_vs_mañana_ex2,"t_m_ex2")

tarde_vs_mañana_ex3 = list("ex3col0t_m", "ex3mojt_m", "ex3velt_m", "ex3bont_m")
venn.cg.enrgo("g1/",tarde_vs_mañana_ex3, "t_m_ex3")

tarde_vs_mañana_ex4 = list("ex4col0t_m", "ex4velt_m")
venn.cg.enrgo("g1/",tarde_vs_mañana_ex4,"t_m_ex4")

#GRUPO 2 ("g2/"):tardes vs mañana entre EXTRACCIONES por cada ecotipo
tarde_vs_mañana_col0 = list("ex1col0t_m", "ex2col0t_m","ex3col0t_m","ex4col0t_m")
venn.cg.enrgo("g2/",tarde_vs_mañana_col0, "t_m_col0")

tarde_vs_mañana_vel = list("ex1velt_m", "ex2velt_m","ex3velt_m","ex4velt_m")
venn.cg.enrgo("g2/",tarde_vs_mañana_vel,"t_m_vel")

tarde_vs_mañana_moj = list("ex1mojt_m", "ex2mojt_m","ex3mojt_m")
venn.cg.enrgo("g2/",tarde_vs_mañana_moj, "t_m_moj")

tarde_vs_mañana_bon = list("ex1bont_m", "ex2bont_m","ex3bont_m")
venn.cg.enrgo("g2/",tarde_vs_mañana_bon,"t_m_bon")

dev.set(dev.prev())
dev.set(dev.next())



#EXTRACCIÓN VS EXTRACCIÓN ANTERIOR:
#GRUPO 3 ("g3/"):EX-EX_anterior cada ecotipo para mañana y tarde
exs_vel_m = list("ex2-ex1velm", "ex3-ex2velm","ex4-ex3velm")
venn.cg.enrgo("g3/",exs_vel_m,"exs_vel_m" )
exs_vel_t = list("ex2-ex1velt", "ex3-ex2velt","ex4-ex3velt")
venn.cg.enrgo("g3/",exs_vel_t,"exs_vel_t")

exs_vel_m1 = list("ex2-ex1velm", "ex3-ex1velm","ex3-ex2velm")
venn.cg.enrgo("g3/",exs_vel_m1, "exs_vel_m1")
exs_vel_t1 = list("ex2-ex1velt", "ex3-ex1velt","ex3-ex2velt")
venn.cg.enrgo("g3/",exs_vel_t1, "exs_vel_t1")

exs_col0_m = list("ex2-ex1col0m", "ex3-ex2col0m","ex4-ex3col0m")
venn.cg.enrgo("g3/",exs_col0_m, "exs_col0_m")
exs_col0_t = list("ex2-ex1col0t", "ex3-ex2col0t","ex4-ex3col0t")
venn.cg.enrgo("g3/",exs_col0_t ,"exs_col0_t")

exs_col0_t1 = list("ex2-ex1col0t", "ex3-ex1col0t","ex3-ex2col0t")
venn.cg.enrgo("g3/",exs_col0_t1 ,"exs_col0_t1")
exs_col0_m1 = list("ex2-ex1col0m", "ex3-ex1col0m","ex3-ex2col0m")
venn.cg.enrgo("g3/",exs_col0_m1 ,"exs_col0_m1")

exs_bon_m = list("ex2-ex1bonm", "ex3-ex1bonm","ex3-ex2bonm") 
venn.cg.enrgo("g3/",exs_bon_m, "exs_bon_m")
exs_bon_t = list("ex2-ex1bont", "ex3-ex1bont","ex3-ex2bont")
venn.cg.enrgo("g3/",exs_bon_t ,"exs_bon_t")

exs_moj_m = list("ex2-ex1mojm", "ex3-ex1mojm","ex3-ex2mojm")
venn.cg.enrgo("g3/",exs_moj_m ,"exs_moj_m")
exs_moj_t = list("ex2-ex1mojt", "ex3-ex1mojt","ex3-ex2mojt")
venn.cg.enrgo("g3/",exs_moj_t ,"exs_moj_t")

#GRUPO 4 ("g4/"):EX-EX_anterior entre ecotipos para mañana y tarde
#MAÑANA
exs2_1_ecot_m = list( "ex2-ex1col0m","ex2-ex1mojm","ex2-ex1velm","ex2-ex1bonm")
venn.cg.enrgo("g4/",exs2_1_ecot_m, "exs2_1_m")

exs3_1_ecot_m = list( "ex3-ex1col0m","ex3-ex1mojm","ex3-ex1velm","ex3-ex1bonm")
venn.cg.enrgo("g4/",exs3_1_ecot_m, "exs3_1_m")

exs3_2_ecot_m = list( "ex3-ex2col0m","ex3-ex2mojm","ex3-ex2velm","ex3-ex2bonm")
venn.cg.enrgo("g4/",exs3_2_ecot_m,"exs3_2_m")

#TARDE
exs2_1_ecot_t = list( "ex2-ex1col0t","ex2-ex1mojt","ex2-ex1velt","ex2-ex1bont")
venn.cg.enrgo("g4/",exs2_1_ecot_t, "exs2_1_t")

exs3_1_ecot_t = list( "ex3-ex1col0t","ex3-ex1mojt","ex3-ex1velt","ex3-ex1bont")
venn.cg.enrgo("g4/",exs3_1_ecot_t,"exs3_1_t")

exs3_2_ecot_t = list( "ex3-ex2col0t","ex3-ex2mojt","ex3-ex2velt","ex3-ex2bont")
venn.cg.enrgo("g4/",exs3_2_ecot_t,"exs3_2_t")


#ECOTIPOS VS MONTEJAQUE:
#GRUPO 5 ("g5/"):Entre las extracciones de cada ecotipo para mañana y tarde:
exs_vel_moj_m = list("ex1vel-mojm", "ex2vel-mojm", "ex3vel-mojm")
venn.cg.enrgo("g5/",exs_vel_moj_m, "vel_moj_m")
exs_vel_moj_t = list("ex1vel-mojt", "ex2vel-mojt", "ex3vel-mojt")
venn.cg.enrgo("g5/",exs_vel_moj_t,"vel_moj_t")

exs_col0_moj_m = list("ex1col0-mojm", "ex2col0-mojm", "ex3col0-mojm")
venn.cg.enrgo("g5/",exs_col0_moj_m,"col0_moj_m")
exs_col0_moj_t = list("ex1col0-mojt", "ex2col0-mojt", "ex3col0-mojt")
venn.cg.enrgo("g5/",exs_col0_moj_t, "col0_moj_t")

exs_bon_moj_m = list("ex1bon-mojm", "ex2bon-mojm", "ex3bon-mojm")
venn.cg.enrgo("g5/",exs_bon_moj_m, "bon_moj_m")
exs_bon_moj_t = list("ex1bon-mojt", "ex2bon-mojt", "ex3bon-mojt")
venn.cg.enrgo("g5/",exs_bon_moj_t, "bon_moj_t")


#GRUPO 6 ("g6/"):Cada extracción entre ecotipos para mañana y tarde
ex1_ecos_moj_m = list("ex1col0-mojm", "ex1vel-mojm","ex1bon-mojm")
venn.cg.enrgo("g6/",ex1_ecos_moj_m, "ecos_moj_ex1m")

ex1_ecos_moj_t = list("ex1col0-mojt", "ex1vel-mojt","ex1bon-mojt")
venn.cg.enrgo("g6/",ex1_ecos_moj_t,"ecos_moj_ex1t")

ex2_ecos_moj_m = list("ex2col0-mojm", "ex2vel-mojm","ex2bon-mojm")
venn.cg.enrgo("g6/",ex2_ecos_moj_m,"ecos_moj_ex2m")

ex2_ecos_moj_t = list("ex2col0-mojt", "ex2vel-mojt","ex2bon-mojt")
venn.cg.enrgo("g6/",ex2_ecos_moj_t,"ecos_moj_ex2t")

ex3_ecos_moj_m = list("ex3col0-mojm", "ex3vel-mojm","ex3bon-mojm")
venn.cg.enrgo("g6/",ex3_ecos_moj_m,"ecos_moj_ex3m")

ex3_ecos_moj_t = list("ex3col0-mojt", "ex3vel-mojt","ex3bon-mojt")
venn.cg.enrgo("g6/",ex3_ecos_moj_t, "ecos_moj_ex3t")


#ECOTIPOS VS COLUMBIA
#GRUPO 7 ("g7/"):Entre las extracciones de cada ecotipo para mañana y tarde:
exs_vel_col0_m = list("ex1vel-col0m", "ex2vel-col0m", "ex3vel-col0m")
venn.cg.enrgo("g7/",exs_vel_col0_m, "vel_col_m")
exs_vel_col0_t = list("ex1vel-col0t", "ex2vel-col0t", "ex3vel-col0t")
venn.cg.enrgo("g7/",exs_vel_col0_t,"vel_col_t")

exs_bon_col0_m = list("ex1bon-col0m", "ex2bon-col0m", "ex3bon-col0m")
venn.cg.enrgo("g7/",exs_bon_col0_m, "bon_col_m")
exs_bon_col0_t = list("ex1bon-col0t", "ex2bon-col0t", "ex3bon-col0t")
venn.cg.enrgo("g7/",exs_bon_col0_t,"bon_col_t")

exs_moj_col0_m = list("ex1moj-col0m", "ex2moj-col0m", "ex3moj-col0m")
venn.cg.enrgo("g7/",exs_moj_col0_m, "moj_col_m")
exs_moj_col0_t = list("ex1moj-col0t", "ex2moj-col0t", "ex3moj-col0t")
venn.cg.enrgo("g7/",exs_moj_col0_t, "moj_col_t")


#GRUPO 8 ("g8/"):Cada extracción entre ecotipos para mañana y tarde
ex1_ecos_col0_m = list("ex1moj-col0m", "ex1vel-col0m","ex1bon-col0m")
venn.cg.enrgo("g8/",ex1_ecos_col0_m, "ecos_col_ex1m")
ex1_ecos_col0_t = list("ex1moj-col0t", "ex1vel-col0t","ex1bon-col0t")
venn.cg.enrgo("g8/",ex1_ecos_col0_t, "ecos_col_ex1t")

ex2_ecos_col0_m = list("ex2moj-col0m", "ex2vel-col0m","ex2bon-col0m")
venn.cg.enrgo("g8/",ex2_ecos_col0_m,"ecos_col_ex2m")
ex2_ecos_col0_t = list("ex2moj-col0t", "ex2vel-col0t","ex2bon-col0t")
venn.cg.enrgo("g8/",ex2_ecos_col0_t, "ecos_col_ex2t")

ex3_ecos_col0_m = list("ex3moj-col0m", "ex3vel-col0m","ex3bon-col0m")
venn.cg.enrgo("g8/",ex3_ecos_col0_m, "ecos_col_ex3m")
ex3_ecos_col0_t = list("ex3moj-col0t", "ex3vel-col0t","ex3bon-col0t")
venn.cg.enrgo("g8/",ex3_ecos_col0_t, "ecos_col_ex3t")



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # UNION DE GENES, SELECCIÓN DE GENES EN MATRIZ DE EXPRESIÓN Y CLUSTERING  # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#Función para seleccionar los genes obtenidos en los contrastes (sin repetirlos)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

lista_genes_tot = function(nom_contr,file_genes_act_tot, file_genes_rep_tot){
  genes_act_tot = vector()
  genes_rep_tot = vector()
  for (i in seq_along(nom_contr)){
    genes_act = read.table((paste(c("results/",nom_contr[i],"_activated",".tsv"),
                                  collapse = "")), as.is = TRUE, header = FALSE)[[1]]
    genes_act_tot = c(union(genes_act, genes_act_tot))
  }
  write.table(genes_act_tot,file =  file_genes_act_tot,
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  
  for (i in seq_along(nom_contr)){
    genes_rep = read.table((paste(c("results/",nom_contr[i],"_repressed",".tsv"),
                                  collapse = "")), as.is = TRUE, header = FALSE)[[1]]
    genes_rep_tot = c(union(genes_rep, genes_rep_tot))
  }
  write.table(genes_rep_tot,file =  file_genes_rep_tot,
              quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  
}

##############################################################################
##############################################################################
################                                            ##################
################               FIN DE LA FUNCIÓN            ##################
################                                            ##################
##############################################################################
##############################################################################



#REDUCCIÓN DE LA MATRIZ DE EXPRESIÓN A LOS GENES DE INTERÉS EN EL ESTUDIO

nomb_genes_act = "results/genes_totales_act.tsv"
nomb_genes_rep = "results/genes_totales_rep.tsv"
nombres_contrastes = list("ex1col0t_m", "ex1mojt_m", "ex1velt_m", "ex1bont_m",
                          "ex2col0t_m", "ex2mojt_m", 
                          "ex2velt_m", "ex2bont_m", "ex3col0t_m", "ex3mojt_m", "ex3velt_m", "ex3bont_m",
                          "ex4col0t_m", "ex4velt_m","ex2-ex1velm", "ex3-ex1velm", "ex3-ex2velm",
                          "ex4-ex1velm","ex4-ex2velm","ex4-ex3velm", "ex2-ex1velt", "ex3-ex1velt", 
                          "ex3-ex2velt","ex4-ex1velt","ex4-ex2velt","ex4-ex3velt", "ex2-ex1col0t", 
                          "ex3-ex1col0t","ex3-ex2col0t","ex4-ex1col0t","ex4-ex2col0t","ex4-ex3col0t",
                          "ex2-ex1col0m", "ex3-ex1col0m","ex3-ex2col0m","ex4-ex1col0m","ex4-ex2col0m",
                          "ex4-ex3col0m","ex2-ex1bonm", "ex3-ex1bonm","ex3-ex2bonm","ex2-ex1mojm", 
                          "ex3-ex1mojm","ex3-ex2mojm", "ex1vel-mojm", "ex2vel-mojm", "ex3vel-mojm",
                          "ex1col0-mojm", "ex2col0-mojm", "ex3col0-mojm","ex1bon-mojm", "ex2bon-mojm", 
                          "ex3bon-mojm","ex1vel-mojt", "ex2vel-mojt", "ex3vel-mojt","ex1col0-mojt", 
                          "ex2col0-mojt", "ex3col0-mojt","ex1bon-mojt", "ex2bon-mojt", "ex3bon-mojt",
                          "ex1vel-col0m", "ex2vel-col0m", "ex3vel-col0m","ex1bon-col0m", "ex2bon-col0m", 
                          "ex3bon-col0m","ex1moj-col0m", "ex2moj-col0m", "ex3moj-col0m","ex1vel-col0t", 
                          "ex2vel-col0t", "ex3vel-col0t", "ex1bon-col0t", "ex2bon-col0t", "ex3bon-col0t",
                          "ex1moj-col0t", "ex2moj-col0t", "ex3moj-col0t")

lista_genes_tot(nombres_contrastes,nomb_genes_act,nomb_genes_rep)

total_act = read.table("results/genes_totales_act.tsv", as.is = TRUE, 
                       header = FALSE)[[1]]
total_rep = read.table("results/genes_totales_rep.tsv", as.is = TRUE, 
                       header = FALSE)[[1]]

summary(total_act)
summary(total_rep)

genes_totales = union(total_act, total_rep)
summary(genes_totales)
head(genes_totales)


#Una vez seleccionados los genes, tenemos que obtener la matriz de expresión para
#estos en lugar de para todos los genes de Arabidopsis.

log.normalized.ecofun.gene.expression <- read.table(file="ecofun/Quantile-normalized.txt",
                                                    header=T,as.is=T)

geneID <- log.normalized.ecofun.gene.expression$geneID
dim(log.normalized.ecofun.gene.expression)
#Quédate sólo con las columnas desde la segunda a la 83 
log.normalized.ecofun.gene.expression <- as.matrix(log.normalized.ecofun.gene.expression[,2:84])
#nombra las filas con los nombres en la columna geneID.
rownames(log.normalized.ecofun.gene.expression) <- geneID
head(log.normalized.ecofun.gene.expression)


reducidos = log.normalized.ecofun.gene.expression[genes_totales,]
dim(reducidos)
head(reducidos)

#aqui es necesario tener cargado el paquete limma
experimental.design = model.matrix(~ -1+factor(c(1,1,1,2,2,2,3,3,3,4,4,4,
                                                 5,5,6,6,6,7,7,7,8,8,8,
                                                 9,9,9,10,10,10,11,11,11,12,12,12,
                                                 13,13,13,14,14,14,15,15,15,16,16,16,
                                                 17,17,17,18,18,18,19,19,19,20,20,20,
                                                 21,21,21,22,22,22,23,23,23,24,24,24,
                                                 25,25,25,26,26,26,27,27,27,28,28,28)))
colnames(experimental.design) = c("ex1velm", "ex1velt", "ex2velm",  "ex2velt",
                                  "ex3velm", "ex3velt",  "ex4velm", "ex4velt",
                                  "ex1bonm", "ex1bont", "ex2bonm", "ex2bont",
                                  "ex3bonm", "ex3bont", "ex1mojm",  "ex1mojt",
                                  "ex2mojm", "ex2mojt",  "ex3mojm",  "ex3mojt",
                                  "ex1col0m","ex1col0t", "ex2col0m", "ex2col0t",
                                  "ex3col0m", "ex3col0t", "ex4col0m", "ex4col0t")


# A continuación, ajustamos la estimación de los niveles de expresión de cada
## gen a un modelo lineal teniendo en cuenta el diseño experimental. Este paso
## fundamentalmente se calcula la media de las réplicas en cada condición.

linear.fit.2 = lmFit(reducidos, experimental.design)
red_gen = linear.fit.2$coefficients
dim(red_gen)
head(red_gen)

#matriz de distancias
library("WGCNA") 
allowWGCNAThreads()

ecofun.correlation = WGCNA::cor(t(x = red_gen))
(dim(x = ecofun.correlation))
ecofun.correlation[1:5, 1:5]

head(ecofun.correlation)

#A partir de la matriz de expresión reducida de estos genes,
#se realizará el clustering


#CLUSTERING AITOE:
library(factoextra) 

k2 <- kmeans((t(red_gen)), centers = 2, nstart = 20)
fviz_cluster(k2, data = (t(red_gen)))

k3 <- kmeans((t(red_gen)), centers = 3, nstart = 20)
fviz_cluster(k3, data = (t(red_gen)))

k4 <- kmeans((t(red_gen)), centers = 4, nstart = 20)
fviz_cluster(k4, data = (t(red_gen)))

#Número óptimo de clusteres
fviz_nbclust(x = red_gen, FUNcluster = kmeans, method = "wss", 
             diss = dist(red_gen, method = "manhattan"))



#a partir de aquí no usarlos.
k5 <- kmeans((t(red_gen)), centers = 5, nstart = 20)
fviz_cluster(k5, data = (t(red_gen)))

k6 <- kmeans((t(red_gen)), centers = 6, nstart = 20)
fviz_cluster(k6, data = (t(red_gen)))

k7 <- kmeans((t(red_gen)), centers = 7, nstart = 20)
fviz_cluster(k7, data = (t(red_gen)))

k8 <- kmeans((t(red_gen)), centers = 8, nstart = 20)
fviz_cluster(k8, data = (t(red_gen)))

k9 <- kmeans((t(red_gen)), centers = 9, nstart = 20)
fviz_cluster(k9, data = (t(red_gen)))

k10 <- kmeans((t(red_gen)), centers = 10, nstart = 20)
fviz_cluster(k10, data = (t(red_gen)))

k11 <- kmeans((t(red_gen)), centers = 11, nstart = 20)
fviz_cluster(k11, data = (t(red_gen)))

k12 <- kmeans((t(red_gen)), centers = 12, nstart = 20)
fviz_cluster(k12, data = (t(red_gen)))


# Hacer clustering solo con los datos de mañana y solo de tarde.
pares = seq(2,28,by=2)
impares = seq(1,27,by=2)

red_gen_m = red_gen[,impares]
red_gen_t = red_gen[,pares]

red_gen_t = t(red_gen_t)
red_gen_t = as.data.frame(red_gen_t)



var0 = apply((t(red_gen_t)), FUN = sd, MARGIN = 2)
no_var = which(var0 == 0)
names = names(no_var)

red_gen_t1 = select(red_gen_t, -AT5G07550)
red_gen_t1 = select(red_gen_t1, -AT1G68875)
red_gen_t1 = select(red_gen_t1, -AT2G19070)

#clustering mañana
k2m <- kmeans((t(red_gen_m)), centers = 2, nstart = 20)
fviz_cluster(k2m, data = (t(red_gen_m)))

k3m <- kmeans((t(red_gen_m)), centers = 3, nstart = 20)
fviz_cluster(k3m, data = (t(red_gen_m)))

k4m <- kmeans((t(red_gen_m)), centers = 4, nstart = 20)
fviz_cluster(k4m, data = (t(red_gen_m)))

#Número óptimo de clusteres
fviz_nbclust(x = red_gen_m, FUNcluster = kmeans, method = "wss", 
             diss = dist(red_gen_m, method = "manhattan"))


#clustering tarde
k2t <- kmeans(red_gen_t1, centers = 2, nstart = 20)
fviz_cluster(k2t, data = red_gen_t1, ellipse = TRUE)

k3t <- kmeans(red_gen_t1, centers = 3, nstart = 20)
fviz_cluster(k3t, data = red_gen_t1, ellipse.type = "confidence")

k4t <- kmeans(red_gen_t1, centers = 4, nstart = 20)
fviz_cluster(k4t, data = red_gen_t1, ellipse.type = "confidence")

#Número óptimo de clusteres
fviz_nbclust(x = red_gen_t1, FUNcluster = kmeans, method = "wss", 
             diss = dist(red_gen_t1, method = "manhattan"))

#Guardo los datos con los que hemos realizado el clustering en ficheros;
write.csv(red_gen,  file = "results/clustering/datos_clustering.csv",row.names = TRUE,
          col.names = TRUE,quote = FALSE)

write.csv(red_gen_m,  file = "results/clustering/datos_clustering_m.csv",row.names = TRUE,
          col.names = TRUE,quote = FALSE)

write.csv(t(red_gen_t1),  file = "results/clustering/datos_clustering_t.csv",row.names = TRUE,
          col.names = TRUE,quote = FALSE)

#heatmap
heatmap(red_gen)








