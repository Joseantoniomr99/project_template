tablaUniprot<- read.csv(file="data.csv", sep="", header = TRUE)
tablaUniprot$Cross.reference<- NULL
tablaUniprot$X.STRING.<-NULL
tablaUniprot[,3] <- gsub(";", tablaUniprot[,3], replacement = "")
tablaproteoma<-  read.table("9606.protein.links.v11.5.txt", sep = "", header = TRUE)
network.table <- openxlsx::read.xlsx("Prey_Lookup_Table.xlsx")
names(network.table)<- c(network.table[1,1],network.table[1,2],network.table[1,3],
                         network.table[1,4],network.table[1,5],network.table[1,6],
                         network.table[1,7],network.table[1,8],network.table[1,9],
                         network.table[1,10],network.table[1,11],network.table[1,12])
network.table<- network.table[-1,]

##Bucle para eliminar los que no tienen valor en codigo ensamble
for(i in 1:length(tablaUniprot$name)){
  if(tablaUniprot$name[i]==""){
    tablaUniprot$name[i]<- NA
  }
}
tablaUniprot<- na.omit(tablaUniprot)

##Buscamos las proteinas con codigo ensamble de la tabla descargadas que aparecen en el proteoma. 
interacc<- tablaUniprot$name %in% tablaproteoma$protein1
v<-c()
for(i in 1:length(interacc)){
  if(interacc[i]==TRUE){
    v<- c(v,i)
  }else{
    i=i+1
  }
}
##Vemos que de 18870 proteinas que tenemos con codigo ensamble encontramos en el proteoma humano 18866, solo perdemos 4 de ellas.
## Una vez que tenemos el numero de indice donde se encuentra pues obtendremos el codigo ensamble.

ens<- c()
uni<- c()
for(i in 1:length(v)){
  ens<-c(ens, tablaUniprot$name[i])
  uni<- c(uni, tablaUniprot$Entry[i])
}
tabla <- data.frame(ens,uni)

##Ahora ya por ultimo toca emplear el mismo metodo para ver cuantos de estos valores los encontramos en nuestra tabla de covid.
matchh<- match(network.table$PreyUniprotAcc,tabla$uni)
PreyUniprotAcc<-c()
ensc<- c()
for (i in 1:length(network.table$PreyUniprotAcc)){
  if((is.na(matchh[i]))==FALSE){
    PreyUniprotAcc<-c(PreyUniprotAcc, network.table$PreyUniprotAcc[i])
  }
}
matchh2<- tabla$uni %in% PreyUniprotAcc

tablaMatch<- data.frame(ensc,PreyUniprotAcc)

for(i in 1:length(matchh2)){
  if(matchh2[i]==TRUE){
    ensc<-c(ensc, tabla$ens[i])
    PreyUniprotAcc<-c(PreyUniprotAcc,tabla$uni[i])
  }
}
tablaMatch<- data.frame(ensc,PreyUniprotAcc)

##Obtenemos la tabla final de genes covid y genes humanos unidos, donde hemos conseguido su codigo ensabmle de string.
## En este caso no perdemos unicamente 4 genes de las proteinas humanas.
## No perdemos nignun gen de covid.  
tablaFinal <- merge(x = network.table, y = tablaMatch, by ="PreyUniprotAcc")
unique(network.table$Bait)
unique(tablaFinal$Bait)

