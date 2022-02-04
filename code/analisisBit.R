library(xlsx)

protTraduc<- read.csv("ProteomaTraducido.csv")
network.table <- openxlsx::read.xlsx("Prey_Lookup_Table.xlsx")
names(network.table)<- c(network.table[1,1],network.table[1,2],network.table[1,3],
                         network.table[1,4],network.table[1,5],network.table[1,6],
                         network.table[1,7],network.table[1,8],network.table[1,9],
                         network.table[1,10],network.table[1,11],network.table[1,12])
network.table<- network.table[-1,]

#Comprobamos cuantos valores de proteinas de nuestra network table se encuentra en nuestro proteoma traducido
# Vemos que solo encontramos 18 de las proteinas, de un total de 332 proteinas. 
ain<- network.table$PreyGeneName %in% protTraduc$protein1
v<-c()
for(i in 1:length(ain)){
  if(ain[i]==TRUE){
    v<- c(v,i)
  }else{
    i=i+1
  }
  
}
 # Añadimos en dos vectores los valores de esas posiciones, y vemos que solo conseguimos 2 genes covid.
prot<- c()
cov<- c()
for(i in 1:length(v)){
  prot<-c(prot, network.table$PreyGeneName[i])
  cov<- c(cov, network.table$Bait[i])
}
