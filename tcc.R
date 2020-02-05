setwd('/home/leonardo/flowmicroarray/imuno/GSE22886_RAW/') # http://bioinformatics.knowledgeblog.org/2011/06/20/analysing-microarray-data-in-bioconductor/

#############
#LOWLEVEL ANALYSIS
#############

library('affy')
library('gcrma')
library('plier')
library('panp')

#datahandling
phenodat = read.csv('phenodata.csv', sep = ',', header = T , row.names = list.celfiles()) #carrega arquivo cel do diretorio
control = data.frame(phenodat[phenodat$treatment == 'Control',], row.names = 2) #manipulacao de dadoshis
celnames = as.character(control$celname) #vetor com nomes dos arquivos cel
CDs = read.csv('CDs.csv', sep = ',', header = T) #lista com IDs dos Clusters Definitions (CDs)
MPs = read.csv('MPs.csv', sep = ',', header = T) #lista com IDs de proteinas de membranas
cytoflow = as.data.frame(read.csv('chartdemo.csv', sep = ',', header = T, row.names = 1, na.strings = T), optional = T) #dados de presenca/ausenca de

#loading cel files
celfiles = ReadAffy(filenames = celnames, phenoData = control, celfile.path = '/home/leonardo/flowmicroarray/imuno/GSE22886_RAW/celfiles/')

#preprocessing
#RMA
celfiles.rma <- rma(celfiles)
eset.rma <- exprs(celfiles.rma)
pa.rma <- pa.calls(celfiles.rma, looseCutoff=0.02, tightCutoff=0.01, verbose = FALSE)
pcalls.rma <- pa.rma$Pcalls
pvalues.rma <- pa.rma$Pvals
##loose cutoff
pa.rma2 <- pa.calls(celfiles.rma, looseCutoff = 0.1, tightCutoff = 0.09, verbose = FALSE)
pcalls.rma2 <- pa.rma2$Pcalls
cds.rma2 = match(as.vector(CDs$probeID), row.names(pcalls.rma2))
cds.rma2 = pcalls.rma2[cds.rma2,]
##cds and membrane protein filtering
cds.rma = match(as.vector(CDs$probeID), row.names(pcalls.rma))
cds.rma = pcalls.gcrma[cds.rma,]
mps.rma = match(as.vector(MPs$probeID), row.names(pcalls.rma))
mps.rma = pcalls.gcrma[mps.rma,]
exp.cdsrma = eset.rma[match(as.vector(CDs$probeID), rownames(eset.rma)),]

###write it
write.table(pcalls.rma, file = 'rmapcalls.csv', sep = ',', col.names = NA)
write.table(pvalues.rma, file = 'rmapvalues.csv', sep = ',', col.names = NA)
write.table(cds.rma, file = 'rmaCDs.csv', sep = ',', col.names = NA)
write.table(mps.rma, file = 'rmaMPs.csv', sep = ',', col.names = NA)


#GCRMA
celfiles.gcrma <- gcrma(celfiles)
eset.gcrma <- exprs(celfiles.gcrma)
pa.gcrma <- pa.calls(celfiles.rma, looseCutoff=0.02, tightCutoff=0.01, verbose = FALSE)
pcalls.gcrma <- pa.gcrma$Pcalls
pvalues.gcrma <- pa.gcrma$Pvals
##cds and membrane protein filtering
cds.gcrma = match(as.vector(CDs$probeID), row.names(pcalls.gcrma))
cds.gcrma = pcalls.gcrma[cds.gcrma,]
mps.gcrma = match(as.vector(MPs$probeID), row.names(pcalls.gcrma))
mps.gcrma = pcalls.gcrma[mps.gcrma,]
exp.cdsgcrma = eset.gcrma[match(as.vector(CDs$probeID), rownames(eset.gcrma)),]
###write it
write.table(pcalls.gcrma, file = 'gcrmapcalls.csv', sep = ',', col.names = NA)
write.table(pvalues.gcrma, file = 'gcrmapvalues.csv', sep = ',', col.names = NA)
write.table(cds.gcrma, file = 'gcrmaCDs.csv', sep = ',', col.names = NA)
write.table(mps.gcrma, file = 'gcrmaMPs.csv', sep = ',', col.names = NA)

#MAS5
celfiles.mas5 <- mas5(celfiles)
eset.mas5 <- exprs(celfiles.mas5)
pa.mas5 <- pa.calls(celfiles.mas5, looseCutoff=0.02, tightCutoff=0.01, verbose = FALSE)
pcalls.mas5 <- pa.mas5$Pcalls
pvalues.mas5 <- pa.mas5$Pvals

##cds filtering
cds.mas5 = match(as.vector(CDs$probeID), row.names(pcalls.mas5))
cds.mas5 = data.frame(pcalls.mas5[cds.mas5,])
mps.mas5 = match(as.vector(MPs$probeID), row.names(pcalls.mas5))
mps.mas5 = pcalls.mas5[mps.mas5,]
exp.cdsmas5 = eset.mas5[match(as.vector(CDs$probeID), rownames(eset.mas5)),]
###write it
write.table(pcalls.mas5, file = 'mas5pcalls.csv', sep = ',', col.names = NA)
write.table(pvalues.mas5, file = 'mas5pvalues.csv', sep = ',', col.names = NA)
write.table(cds.mas5, file = 'mas5CDs.csv', sep = ',', col.names = NA)
write.table(mps.mas5, file = 'mas5MPs.csv', sep = ',', col.names = NA)

#PLIER
celfiles.plier <- justPlier(celfiles, normalize = TRUE)
eset.plier <- exprs(celfiles.plier)
pa.plier <- pa.calls(celfiles.plier, looseCutoff=0.02, tightCutoff=0.01, verbose = FALSE)
pcalls.plier <- pa.plier$Pcalls
pvalues.plier <- pa.plier$Pvals
##cds filtering
cds.plier = match(as.vector(CDs$probeID), row.names(pcalls.plier))
cds.plier = pcalls.plier[cds.plier,]
mps.plier = match(as.vector(MPs$probeID), row.names(pcalls.plier))
mps.plier = pcalls.plier[mps.plier,]
exp.cdsplier = eset.plier[match(as.vector(CDs$probeID), rownames(eset.plier)),]
###write it
write.table(pcalls.plier, file = 'mas5pcalls.csv', sep = ',', col.names = NA)
write.table(pvalues.plier, file = 'mas5pvalues.csv', sep = ',', col.names = NA)
write.table(cds.plier, file = 'mas5CDs.csv', sep = ',', col.names = NA)
write.table(mps.plier, file = 'mas5MPs.csv', sep = ',', col.names = NA)

######### Controle de Qualidade
library(affyPLM)

par(mfrow = c(2,2))
par(mar = c(7,4,4,2)+0.1)
boxplot(celfiles, col = 'grey50', ylab = 'Sinal de Intensidade das Sondas - (Log2)', las = 2, names = phenodat$cell)
boxplot(log2(eset.mas5), col = 'grey50', main= "MAS5", ylab = 'Sinal de Intensidade das Sondas - (Log2)', xlab = '', las = 2, names = phenodat$cell, ylim= c(-5, 15))
boxplot(eset.plier, col = 'grey50', main= "PLIER", xlab = '',las = 2, names = phenodat$cell, ylim= c(-5, 15))
boxplot(eset.rma, col = 'grey50', main = 'RMA', ylab = 'Sinal de Intensidade das Sondas - (Log2)', las = 2, names = phenodat$cell, ylim= c(2, 16))
boxplot(eset.gcrma, col = 'grey50', main= "GC-RMA", xlab = '',las = 2, names = phenodat$cell, ylim= c(2, 16))
dev.off()

celfiles.qc <- fitPLM(celfiles)

par(mfrow = c(1,2))
par(mar = c(7,4,4,2)+0.1)
RLE(celfiles.qc, main="RLE", las = 2, names = phenodat$cell)
NUSE(celfiles.qc, main="NUSE", las = 2, names = phenodat$cell)
dev.off()

#############
#HIGH LEVEL ANALYSIS
#############

typecell = function(j){ #extrai os nomes das amostras por categorias celular
  as.vector(rownames(control)[control$type == colnames(cytoflow)[j]])
}

extractvexp = function(a,j, method){ #extrai os valores de expressao por categoria celular
  if(method == 1){log2(exp.cdsmas5[which(cytoflow[j] == a), typecell(j)])}
  else {if(method == 2){exp.cdsplier[(cytoflow[j] == a), typecell(j)]}
    else{if(method == 3){exp.cdsrma[(cytoflow[j] == a), typecell(j)]}
      else{if(method == 4){exp.cdsgcrma[(cytoflow[j] == a), typecell(j)]}}}}
}

extractpanpcall = function(i,j, method){ #extrai as absoluts calls por categoria celular, baseado em PANP
  if(method == 1){cds.mas5[i, typecell(j)]}
  else {if(method == 2){cds.plier[i, typecell(j)]}
    else{if(method == 3){cds.rma[i, typecell(j)]}
      else{if(method == 4){cds.gcrma[i, typecell(j)]}}}}
}

par(mfrow = c(2,2))
boxplot(extractvexp('P', 1, 1)[,4], col = 'grey50', extractvexp('A', 1, 1)[,4], main = 'MAS5', names = c('P','A'), ylab= "Nivel de Expressao - (Log2)", xlab = 'Cluster Definition Calls (CDs)', ylim = c(-5,14))
abline(h=bmas)
boxplot(extractvexp('P', 1, 2)[,4], col = 'grey50', extractvexp('A', 1, 2)[,4], main = 'PLIER', names = c('P','A'), xlab = 'Cluster Definition Calls (CDs)', ylim = c(-5,14))
abline(h=bplier)
boxplot(extractvexp('P', 1, 3)[,4], col = 'grey50',extractvexp('A', 1, 3)[,4], main = 'RMA', names = c('P','A'), ylab= "Nivel de Expressao - (Log2)", xlab = 'Cluster Definition Calls (CDs)', ylim = c(2,14))
abline(h=brma)
boxplot(extractvexp('P', 1, 4)[,4], col = 'grey50', extractvexp('A', 1, 4)[,4], main = 'GCRMA', names = c('P','A'), xlab = 'Cluster Definition Calls (CDs)', ylim = c(2,14))
abline(h=bgcrma)
dev.off()

tcc1 <- function(j,s,method){
  as.data.frame(ifelse(cytoflow[,j] == 'NA' & extractpanpcall(,j,method)[,s] == 'P', 'UP', 
         ifelse(cytoflow[,j] == 'NA' & extractpanpcall(,j,method)[,s] == 'A', 'UA', 
                ifelse(cytoflow[,j] == 'P' & extractpanpcall(,j,method)[,s] == 'P', 'PP', 
                       ifelse(cytoflow[,j] == 'A' & extractpanpcall(,j,method)[,s] == 'A', 'AA', 
                              ifelse(cytoflow[,j] == 'A' & extractpanpcall(,j,method)[,s] == 'P', 'AP', 
                                     ifelse(cytoflow[,j] == 'P' & extractpanpcall(,j,method)[,s] == 'A', 'PA', 'NA')))))))
}

tcc2 <- function(j,b,s,method){
  as.data.frame(ifelse(cytoflow[,j] == 'NA' & makecall(b,s,method) == 'P', 'UP', 
         ifelse(cytoflow[,j] == 'NA' & makecall(b,s,method) == 'A', 'UA', 
                ifelse(cytoflow[,j] == 'P' & makecall(b,s,j) == 'P', 'PP', 
                       ifelse(cytoflow[,j] == 'A' & makecall(b,s,method) == 'A', 'AA', 
                              ifelse(cytoflow[,j] == 'A' & makecall(b,s,method) == 'P', 'AP', 
                                     ifelse(cytoflow[,j] == 'P' & makecall(b,s,method) == 'A', 'PA', 'NA')))))))
}

summary(extractpanpcall(,1,1)[4] == as.data.frame(cytoflow[,1]))
summary(tcc2(1,bmas,4))
dim(makecall(bmas,1))
##########HIGHER LEVEL ANALYSIS

###PLOTS

par(mfrow= c(2,2),mar = c(7,4,4,2) + 0.1)
boxplot(extractvexp('P',1,1), main = 'Cluster Definition Call P', col = 'grey50', ylab = 'Sinal de Intensidade - (MAS5)', las = 2, names = control[typecell(1),3], ylim= c(0,14))
abline(h=bmas)
boxplot(extractvexp('A',1,1), main = 'Cluster Definition Call A', col = 'grey50', las = 2, names = control[typecell(1),3], ylim= c(0,14))
abline(h=bmas)

boxplot(extractvexp('P',1,2), main = 'Cluster Definition Call P', col = 'grey50', ylab = 'Sinal de Intensidade - (PLIER)', las = 2, names = control[typecell(1),3], ylim= c(0,14))
abline(h=bplier)
boxplot(extractvexp('A',1,2), main = 'Cluster Definition Call A', col = 'grey50', las = 2, names = control[typecell(1),3], ylim= c(0,14))
abline(h=bplier)
dev.off()

par(mfrow= c(2,2),mar = c(7,4,4,2) + 0.1)
boxplot(extractvexp('P',1,3), main = 'Cluster Definition Call P', col = 'grey50', ylab = 'Sinal de Intensidade - (RMA)', las = 2, names = control[typecell(1),3], ylim= c(0,14))
abline(h=brma)
boxplot(extractvexp('A',1,3), main = 'Cluster Definition Call A', col = 'grey50', las = 2, names = control[typecell(1),3], ylim= c(0,14))
abline(h=brma)

boxplot(extractvexp('P',1,4), main = 'Cluster Definition Call P', col = 'grey50', ylab = 'Sinal de Intensidade - (GCRMA)', las = 2, names = control[typecell(1),3], ylim= c(0,14))
abline(h=bgcrma)
boxplot(extractvexp('A',1,4), main = 'Cluster Definition Call A', col = 'grey50', las = 2, names = control[typecell(1),3], ylim= c(0,14))
abline(h=bgcrma)
dev.off()
#Gerando o cutoff e determinando detection call

makecall = function(b,s,method){ #funcao alternativa para definicao de absolut call, baseado em valor escolhido
  if(method == 1){ifelse(log2(exp.cdsmas5[,s]) >= b, 'P', 'A')}
  else {if(method == 2){ifelse(exp.cdsplier[,s] >= b, 'P', 'A')}
    else{if(method == 3){ifelse(exp.cdsrma[,s] >= b, 'P', 'A')}
      else{if(method == 4){ifelse(exp.cdsgcrma[,s] >= b, 'P', 'A')}}}}
}

#todos para tcell 4amostra
bmas <- ((quantile(extractvexp('P',1,1))[2] + quantile(extractvexp('A',1,1))[4])/2)
dcallmas <- as.data.frame(makecall(bmas,1)[, typecell(1)[4]])
Pmas <- rownames(dcallmas)[dcallmas == 'P']
Amas <- rownames(dcallmas)[dcallmas == 'A']

bplier <- (quantile(extractvexp('P',1,2))[2] + quantile(extractvexp('A',1,2))[4])/2
dcallplier <- as.data.frame(makecall(bplier,2)[, typecell(1)[4]])
Pplier <- rownames(dcallplier)[dcallplier == 'P']
Aplier <- rownames(dcallplier)[dcallplier == 'A']


brma <- (quantile(extractvexp('P',1,3))[2] + quantile(extractvexp('A',1,3))[4])/2
dcallrma <- as.data.frame(makecall(brma,3)[, typecell(1)[4]])
Prma <- rownames(dcallrma)[dcallrma == 'P']
Arma <- rownames(dcallrma)[dcallrma == 'A']

bgcrma <- (quantile(extractvexp('P',1,4))[2] + quantile(extractvexp('A',1,4))[4])/2
dcallgcrma <- as.data.frame(makecall(bgcrma,4)[, typecell(1)[4]])
Pgcrma <- rownames(dcallgcrma)[dcallgcrma == 'P']
Agcrma <- rownames(dcallgcrma)[dcallgcrma == 'A']

Pt <- as.data.frame(as.vector(c('MAS5',Pmas,'PLIER', Pplier,'RMA', Prma, 'GCRMA', Pgcrma)))
At <- as.data.frame(as.vector(c('MAS5',Amas,'PLIER', Aplier,'RMA', Arma, 'GCRMA', Agcrma)))

write.csv(Pt, file = 'Pt.csv')
write.csv(At, file = 'At.csv')

cd8a <- as.data.frame(tcc1(1,1)[,5])
cd8b <- as.data.frame(tcc2(1,bmas,5))

cd4a <- as.data.frame(tcc1(1,1)[,5])
cd4b <- as.data.frame(tcc2(1,bmas,5))

##remendos

as.data.frame(rownames(cd8a)[which(cd8a == 'UP')])
as.data.frame(rownames(cd8b)[which(cd8b == 'UP')])
as.data.frame(rownames(cd4a)[which(cd4a == 'UP')])
as.data.frame(rownames(cd4b)[which(cd4b == 'UP')])

a=rownames(cd8a)[which(cd8a == 'UP')]
b=rownames(cd8b)[which(cd8b == 'UP')]
c=rownames(cd4a)[which(cd4a == 'UP')]
d=rownames(cd4b)[which(cd4b == 'UP')]

write.csv(c('',a,'',b,'',c,'',d), file = 'up.csv')

e=rownames(cd8a)[which(cd8a == 'UA')]
f=rownames(cd8b)[which(cd8b == 'UA')]
g=rownames(cd4a)[which(cd4a == 'UA')]
h=rownames(cd4b)[which(cd4b == 'UA')]


write.csv(c('',e,'',f,'',g,'',h), file = "ua.csv")


summary(tcc1(j = 1, s = 4, method = 1))
summary(tcc2(j = 1, b = bmas, s=4, method = 1 ))
summary(tcc1(j = 1, s = 4, method = 2))
summary(tcc2(j = 1, b = bplier, s=4, method = 2 ))
summary(tcc1(j = 1, s = 4, method = 3))
summary(tcc2(j = 1, b = brma, s=4, method = 3 ))
summary(tcc1(j = 1, s = 4, method = 4))
summary(tcc2(j = 1, b = bgcrma, s=4, method = 4 ))

summary(tcc1(j = 1, s = 5, method = 1))
summary(tcc2(j = 1, b = bmas, s=4, method = 1 ))
summary(tcc1(j = 1, s = 5, method = 2))
summary(tcc2(j = 1, b = bplier, s=4, method = 2 ))
summary(tcc1(j = 1, s = 5, method = 3))
summary(tcc2(j = 1, b = brma, s=4, method = 3 ))
summary(tcc1(j = 1, s = 5, method = 4))
summary(tcc2(j = 1, b = bgcrma, s=4, method = 4 ))

summary(tcc2(j = 1, b = bgcrma, s=4, method = 4 ))
summary(tcc2(j = 1, b = bgcrma, s=5, method = 4 ))

cd8gc <- tcc2(j = 1, b = bgcrma, s=4, method = 4 )
cd4gc <- tcc2(j = 1, b = bgcrma, s=5, method = 4 )
a1<- rownames(cd8gc)[which(cd8gc == 'UP')]
a2<- rownames(cd4gc)[which(cd4gc == 'UP')]
b1<- rownames(cd8gc)[which(cd8gc == 'UA')]
b2<- rownames(cd4gc)[which(cd4gc == 'UA')]
write.csv(c('',a1,'',a2,'',b1,'',b2), 'gccut.csv')
cd8gc2 <- tcc1(j = 1, s=4, method = 4 )
cd4gc2 <- tcc1(j = 1, s=5, method = 4 )
a1<- rownames(cd8gc2)[which(cd8gc2 == 'UP')]
a2<- rownames(cd4gc2)[which(cd4gc2 == 'UP')]
b1<- rownames(cd8gc2)[which(cd8gc2 == 'UA')]
b2<- rownames(cd4gc2)[which(cd4gc2 == 'UA')]
write.csv(c('',a1,'',a2,'',b1,'',b2), 'gccut2.csv')

a1cut <- rownames(cd8gc)[which(cd8gc == 'UP')]
a1panp <- rownames(cd8gc2)[which(cd8gc2 == 'UP')]
b1cut <- rownames(cd8gc)[which(cd8gc == 'UA')]
b1panp <- rownames(cd8gc2)[which(cd8gc2 == 'UA')]
write.csv(c('',a1cut,'',a1panp,'',b1cut,'',b1panp), 'gccutpanp.csv')
