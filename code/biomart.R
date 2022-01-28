library(biomaRt)
library(xlsx)
library(dplyr)
library(gtools)
##Leemos la tabla de estudio y hacemos los cambios necesarios ya que no se hace una lectura correcta.
network.table <- openxlsx::read.xlsx("Prey_Lookup_Table.xlsx")
names(network.table)<- c(network.table[1,1],network.table[1,2],network.table[1,3],
                         network.table[1,4],network.table[1,5],network.table[1,6],
                         network.table[1,7],network.table[1,8],network.table[1,9],
                         network.table[1,10],network.table[1,11],network.table[1,12])
network.table<- network.table[-1,]

##Cambiamos de uniprot a ensamble usando el paquete biomart. Tomaremos la columna con código uniprot de nuestra tabla de estudio y la cambiaremos a ensamble, que es el formato que tiene la tabla del proteoma.
ensembl <- useEnsembl(biomart = "genes")
ensembl <- useDataset(dataset = "hsapiens_gene_ensembl", mart = ensembl)
listAttributes(ensembl)[grep("uniprot", listAttributes(ensembl)[,1]), ]
uniprot_vector<- network.table$PreyUniprotAcc
codigoUniProt<- getBM(attributes=c("uniprotswissprot", "ensembl_peptide_id"),filter="uniprotswissprot", 
      values=uniprot_vector, mart=ensembl)

##Leemos el proteoma humano y le quitamos los numeros de delante, para asi quedarnos unicamente con el codigo ensamble.
proteoma.table <- read.table("9606.protein.links.v11.5.txt", sep = "", header = TRUE)
proteoma.table[,1] <- gsub("9606.", proteoma.table[,1], replacement = "")
proteoma.table[,2] <- gsub("9606.", proteoma.table[,2], replacement = "")

##Guardaremos en nuestra variable ain true o false según se encuentre ese código ensamble en el proteoma o no.
##Guardaremos en el vector v, las distintas posiciones donde haya un valor de true, es decir se encuentra el código ensamble.
ain<- codigoUniProt$ensembl_peptide_id %in% proteoma.table$protein1
v<-c()
for(i in 1:length(ain)){
  if(ain[i]==TRUE){
    v<- c(v,i)
  }else{
    i=i+1
  }
    
}

##Crearemos dos vectores un vector ens y otro uni, que almacenaran el valor uniprot y ensable en aquellas posiciones cuyo match con el proteoma era true.
ens<- c()
uni<- c()
for(i in 1:length(v)){
      ens<-c(ens, codigoUniProt$ensembl_peptide_id[i])
      uni<- c(uni, codigoUniProt$uniprotswissprot[i])
}

##Comprobamos que todos nuestros codigos uniprot obtenidos por biomart se encuentran en mi tabla dada.
codigoUniProt$uniprotswissprot%in%network.table$PreyUniprotAcc

##Todos se han conseguido traducir a ensamble, pero el problema es que no todos los ensambles están en el proteoma.
##Solo los conseguidos en el vector ens. Veremos a que codigo de uniprot pertenecen. Crearemos un dataframe.
tablaCambio <- data.frame(ens,uni)

#PreyUniprotAcc<- unique(PreyUniprotAcc)
#unico<- unique(network.table$PreyUniprotAcc)


##DE 332 proteinas diferentes solo conseguimos el codigo ensamble de 107 de ellas.
##Vamos a quedarnos solo con uno de los codigos ensamble que tenemos.
PreyUniprotAcc<- c()
diferEnsamb<- c()
for(i in 1:length(tablaCambio$uni)){
    if( (tablaCambio$uni[i] %in% PreyUniprotAcc)==FALSE ){
      PreyUniprotAcc <- c(PreyUniprotAcc, tablaCambio$uni[i])
       diferEnsamb <- c(diferEnsamb, tablaCambio$ens[i])
  }
}

##Unimos los dos dataframe
tablaCambio<- data.frame(diferEnsamb,PreyUniprotAcc)
final <- merge(x = tablaCambio, y = network.table, by ="PreyUniprotAcc")


##Veremos cuantos genes de covid perdemos. Perdemos 4 genes de covid.
## Ya que perdemos 225 de los genes humanos, pues alguno de estos están unidos a esos 4 de covid. 
## Y por consecuencia también se pierden 
unique(network.table$Bait)
unique(final$Bait)
