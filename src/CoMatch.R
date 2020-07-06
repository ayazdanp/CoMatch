############
##PLS 
##Aida Yazdanparast
##12/2/2018
############

##function to prepare input data for PLS for all drugs    

##read data

setwd("H:/research/multi-omics/PLS model")

setwd("/Users/ayazdanpar/Downloads/multiomics_paper/")


if (!require(xlsx)){install.packages("xlsx")}

#tissue data survical result
 tumor_surv = read.xlsx("data/tumor_surv_result.xlsx" , sheetIndex =1)

cell_drug = read.csv("data/cell_line_drug_result.csv")

#save.image("Pls.RData")

load("Pls.RData")

##function to read all cell sreadsheets into list
readfunc=function(data,path){
  list.data=list()
  for(i in 1:5)
    list.data[[i]] = read.table(paste("data/cell-drug result/",data, i,path , sep=""), sep="\t", header=T )
  return(list.data)
}

list.GE=list()
list.GE= readfunc(data="gene",path="-cell-drug prediction.txt")

list.GE = sapply(list.GE, function(x) {row.names(x) <- x[,1] ; x[,-1]} , simplify = F )

list.CNV=readfunc(data="CNV", path="-cell-drug prediction.txt")

readfunc1=function(data){
  library(xlsx)
  list.data=list()
  for(i in 1:5)
    list.data[[i]] = read.xlsx(paste("data/cell-drug result/",data , sep=""), sheetIndex=i )
  return(list.data)
}

list.mut=readfunc1(data="mutation cell-drug prediction.xlsx")

##function to merge each level data with tissue
# 
# function(tumoredata, list.data){
#    sapply(list.data, function(x) merge(x, tumoredata, by.y="Row.names" , by.x= ) )   
# }


##For gene cell line data: pick only probes with lowest p-values
#this is with only cell line data

##1. read data for all drugs
##2. exclude probes
##3. merge with tissue

#gene data
ge.data = read.table("data/breast_ge.tsv", sep=" " , header=T, skip = 1)
row.names(ge.data) = paste(ge.data$X54675, ge.data$X54675.1 , sep="_")

uniqfunc = function(x){ 
  
  t = merge(ge.data[,1:2], x , by="row.names")  
  rownames(t) = t[,1]
  t = t [,-1]
  genes.cell = unique(t[,1])       
  
  d = list()
  s=list()
  library("dplyr") 
  for(j in 1:length(genes.cell)){
    d[[j]] =filter(t, t[,1] %in% genes.cell[j])
    s[[j]] = filter(d[[j]], d[[j]][,4] == min(unlist(d[[j]][,4])))
  }
  return(s)
}

ge.unique=list()

ge.unique= sapply(list.GE, uniqfunc, simplify = F)

#data of selected 

cell_drug_unique = sapply(ge.unique, function(x) bind_rows(x), simplify = F)

cell_drug_unique = sapply(cell_drug_unique, function(x) { row.names(x) <- paste(x$X54675 , x$X54675.1 , sep="_"); x}, simplify = F)


#merging gene, CNV and mutation with tissue

geneTC_unique =sapply(cell_drug_unique, function(x) merge(x[,-2], tumor_surv[,c(1,4,3)], by.x="X54675", by.y="Row.names"),simplify=F)


Colnames = c("gene", "cell.beta" , "cell.P", "tissue.beta",
             "tissue.P")

geneTC_unique= lapply(geneTC_unique,setNames, Colnames)

##merge mutation and cnv list with tissue then merge them with gene list

cnvTC_unique=sapply(list.CNV, function(x) merge(x, tumor_surv[,c(1,7,6)], by.x="X", by.y="Row.names"),simplify=F)

mutTC_unique=sapply(list.mut, function(x) merge(x, tumor_surv[,c(1,9,8)], by.x="NA.", by.y="Row.names"),simplify=F)

cnvTC_unique= lapply(cnvTC_unique,setNames, Colnames)
mutTC_unique= lapply(mutTC_unique,setNames, Colnames)

##take significant p-values

geneTC_unique_sig=sapply(geneTC_unique,function(x) {x[which(x$cell.P >0.05 |x$tissue.P >0.05),c(2,3,4,5)] <- NA; x}, simplify=F)
cnvTC_unique_sig=sapply(cnvTC_unique,function(x) {x[which(x$cell.P >0.05 |x$tissue.P >0.05),c(2,3,4,5)] <- NA; x}, simplify=F)
mutTC_unique_sig=sapply(mutTC_unique,function(x) {x[which(x$cell.P >0.05 |x$tissue.P >0.05),c(2,3,4,5)] <- NA; x}, simplify=F)

geneTC_unique_sig = sapply(geneTC_unique_sig , function(x) {x$gene = paste(x$gene , "GE" , sep=".") ; x}, simplify = F)
cnvTC_unique_sig = sapply(cnvTC_unique_sig , function(x) {x$gene = paste(x$gene , "cnv" , sep=".") ; x}, simplify = F)
mutTC_unique_sig = sapply(mutTC_unique_sig , function(x) {x$gene = paste(x$gene , "mut" , sep=".") ; x}, simplify = F)

##mergind all data sets together

coefficients= list()
for(i in 1:5){
  coefficients[[i]] = do.call("rbind", list(geneTC_unique_sig[[i]],cnvTC_unique_sig[[i]],mutTC_unique_sig[[i]]))
}

coefficients.comp = sapply(coefficients, function(x) {x=x[complete.cases(x),] ; x}, simplify = F)
coefficients.comp = lapply(coefficients.comp,  function(x) {x$sig.cell.beta= sign(x$cell.beta) ;
x$sig.tissue.beta = sign(x$tissue.beta);
x})

##PLS modeling

library(pls)
pls.result <- sapply(coefficients.comp, function(x) plsr( c(tissue.beta*c(1-tissue.P)) ~ c(cell.beta*(1-cell.P))
                                                          , data = x ,validation = "LOO" , scale = TRUE),
                     simplify = F)

#cor(pls.result[[5]]$scores), abs(pls.result[[5]]$Yscores))

scores.pls  =sapply(pls.result, function(x) {t= cbind(x$scores , x$Yscores); t} , simplify = F)


# #spls fit
# library("spls")
#   
# cv_input1 <- cv.spls( x=c(test2$cell.beta*(1-test2$cell.P)), y=c(test2$tissue.beta*c(1-test2$tissue.P)), 
#                       eta = seq(0.1,0.9,0.1), K = c(1) )
# 
# f_input1 <- spls( x=c(test2$cell.beta*(1-test2$cell.P)), y=c(test2$tissue.beta*c(1-test2$tissue.P)), 
#                   K=1, eta=cv_input1$eta.opt , 
#                   scale.x=T,
#                   scale.y = T)
# 

five.scores = list()

five.scores$gemcitabine.scores = as.data.frame(cbind(as.character(coefficients.comp[[1]][,1]), scores.pls[[1]]))
five.scores$doxorubicin.scores = as.data.frame(cbind(as.character(coefficients.comp[[2]][,1]), scores.pls[[2]]))
five.scores$cyclophosphamide.scores = as.data.frame(cbind(as.character(coefficients.comp[[3]][,1]), scores.pls[[3]]))
five.scores$paclitaxel.scores = as.data.frame(cbind(as.character(coefficients.comp[[4]][,1]), scores.pls[[4]]))
five.scores$docetaxel.scores = as.data.frame(cbind(as.character(coefficients.comp[[5]][,1]), scores.pls[[5]]))

###Tamoxifen , doxorubicin

five.scores = sapply(five.scores , function(x) { colnames(x) = c("gene" , "x.score" , "y.score") ; x} , simplify = F)
five.scores = sapply(five.scores, function(x) { x$x.score =abs(as.numeric(as.character(x$x.score))); 
x$y.score =abs(as.numeric(as.character(x$y.score))) ; x }
, simplify = F)

plot(abs(five.scores$gemcitabine.scores$x.score) , abs(five.scores$gemcitabine.scores$y.score))

# five.scores$gemcitabine.scores[order(abs(five.scores$gemcitabine.scores$y.score), decreasing = T), ]
# five.scores$gemcitabine.scores[order(abs(five.scores$gemcitabine.scores$x.score), decreasing = T), ]
# 

##looking a some genes in doxorubicin, these are the genes with high scores in both x and y 
##adding the beta signs to result

five.scores1=list()
five.scores1$gemcitabine.scores= merge(five.scores[[1]], coefficients.comp[[1]][,c(1, 6,7)] , by="gene")
five.scores1$doxorubicin.scores= merge(five.scores[[2]], coefficients.comp[[2]][,c(1, 6,7)] , by="gene" )    
five.scores1$cyclophosphamide.scores = merge(five.scores[[3]], coefficients.comp[[3]][,c(1, 6,7)] , by="gene" )    
five.scores1$paclitaxel.scores  = merge(five.scores[[4]], coefficients.comp[[4]][,c(1, 6,7)] , by="gene" )    
five.scores1$docetaxel.scores = merge(five.scores[[5]], coefficients.comp[[5]][,c(1, 6,7)] , by="gene" )    


write.xlsx(  five.scores$doxorubicin.scores , "result/five.scores.xlsx", sheetName = "doxorubicin.scores",append = T)

writeFun= function(name){
  t =   write.xlsx(  five.scores1$name , "result/five.scores1.xlsx", sheetName = paste(name),append = T)
  return(t)
}

writeFun(gemcitabine.scores)

write.xlsx(  five.scores1$gemcitabine.scores , "result/five.scores1.xlsx", sheetName = "gemcitabine.scores",append = T)
write.xlsx(  five.scores1$doxorubicin.scores , "result/five.scores1.xlsx", sheetName = "doxorubicin.scores",append = T)
write.xlsx(  five.scores1$cyclophosphamide.scores , "result/five.scores1.xlsx", sheetName = "cyclophosphamide.scores",append = T)
write.xlsx(  five.scores1$paclitaxel.scores , "result/five.scores1.xlsx", sheetName = "paclitaxel.scores",append = T)
write.xlsx(  five.scores1$docetaxel.scores , "result/five.scores1.xlsx", sheetName = "docetaxel.scores",append = T)

###########
###plot the vectors


library(pls)
pls.result <- sapply(coefficients.comp, function(x) plsr( c(tissue.beta*c(1-tissue.P)) ~ c(cell.beta*(1-cell.P))
                                                          , data = x ,validation = "LOO" , scale = TRUE),
                     simplify = F)

doxo.coeff=coefficients.comp[[2]]

which(doxo.coeff$cell.beta* (1-doxo.coeff$cell.P)<(-3) & doxo.coeff$tissue.beta*c(1-doxo.coeff$tissue.P) >0 )


labels.gene = c("ABCB11.GE",
                "C15orf54.GE",
                "BSND.GE",
                "CCDC33.GE",
                "GABRA5.GE",
                "SLC6A2.GE",
                "IL3RA.GE",
                "HERC2.GE",
                "POPDC2.GE",
                "FSTL4.GE",
                "SEMA4F.GE",
                "MYO3A.mut",
                "MTOR.cnv",
                "ANGPTL7.cnv",
                "FBXO2.cnv",
                "AGTRAP.cnv",
                "MTHFR.cnv",
                "MAD2L2.cnv",
                "FBXO6.cnv",
                "FBXO44.cnv",
                "PTCHD2.cnv",
                "CACNA1B.GE",
                "NPAT.GE",
                "PTGDR.GE")

row.names(doxo.coeff)=doxo.coeff$gene
partial=doxo.coeff[row.names(doxo.coeff)%in% labels.gene , ]

ggplot(doxo.coeff)+
  geom_point(aes(x=c(cell.beta* (1-cell.P)), y=(tissue.beta*c(1-tissue.P)) ,colour = row.names(doxo.coeff) %in% labels.gene))+
  scale_colour_manual(values = c("blue", "red"))+
  geom_hline(yintercept=0, color="black")+
  geom_vline(xintercept = 0 , color="black")+
  geom_text(data=partial,
            aes(x=c(cell.beta* (1-cell.P)), y=(tissue.beta*c(1-tissue.P)),label=row.names(partial)))+
  theme_bw()+
  xlab("Cell line Coeff vector")+
  ylab("Tissue Coeff vector")+ theme(axis.line = element_line(colour = "black"),
                                     panel.grid.major = element_blank(),
                                     panel.grid.minor = element_blank(),
                                     panel.border = element_blank(),
                                     panel.background = element_blank()) 



geom_text(aes(x=c(cell.beta* (1-cell.P)), y=(tissue.beta*c(1-tissue.P)) , label=row.names(doxo.coeff)))

geom_text(aes(x=c(cell.beta* (1-cell.P)), y=(tissue.beta*c(1-tissue.P)) , label=row.names(doxo.coeff)[row.names(doxo.coeff)%in% labels.gene]))

