
###### PASOS DEL ALGORITMO
#Matriz adyacencia donde los 0 cambian a 1 cuando interacciona. Poner que el grafo no sea dirigido y evitamos bidireccionalidad.
#0. Definición de los datos
library(xlsx)
library(igraph)
library(dplyr)
library(gtools)
prot1<- c("M6P3","RALA","GABARAPL2","GABARAPL2","GABARAPL2","GABARAPL2","GABARAPL2","OTC", "FUNDC1","PLIN3", "FUNDC1")
prot2<- c("PLIN3","EXOC8","CD300C","UVRAG","FUNDC1","PIK3R4","ATG14","ARG2","OTC", "RALA", "EXOC8")
covid1<-c("E","E","E","M","M")
protcov<-c("M6P3","EXOC8","GABARAPL2","PIK3R4","ARG2")
humana<- data.frame(prot1,prot2)
covid<- data.frame(covid1,protcov)
humana.igraph <- graph_from_data_frame(humana)
covid.igraph <- graph_from_data_frame(covid, directed = FALSE)

#1. Obtener combinaciones posibles de las proteinas unidas al covid
covprot<- covid$protcov
N<- length(covprot)
n<-2
c<- combinations(N,n,covprot)

#2. Coger estas combinaciones y buscar las rutas entre ellas usando el proteoma humano.
rutas<- c()
for(i in 1:length(c[,1])){
  for(j in 1:length(c[,2])){
    camino<- all_simple_paths(humana.igraph, from=c[i,1], to=c[j,2], mode = "all") 
    rutas <- c(rutas ,camino)
  }
}
#c<- c()
#for( i in 1:length(lista)){
#c<- c(make_graph(lista[[i]],c))}


#3. Coger y formar combinaciones de nodos 2 a 2. Tomar la conexion entre una proteina y la adjunta.

vS1<- c()
vS2<- c()
for(i in 1:length(rutas)){
  for(j in 2:length(rutas[[i]])-1){
    ruta1<- rutas[[i]][j]$name
    ruta2<- rutas[[i]][j+1]$name
    vS1<- c(vS1, ruta1)
    vS2 <- c(vS2,ruta2 )
  }
}

# Guardamos combinaciones en el dataframe, para obtener la tabla de relaciones.

rutaFinal<- data.frame(vS1,vS2)

#4.Borramos los datos repetidos, ya que habran filas que se repiten.
rutaFinal<- rutaFinal%>%distinct()
rutaFinalgrafo<- graph_from_data_frame(rutaFinal, directed = FALSE)
plot(rutaFinalgrafo)

#5. Unimos con proteinas covid, podemos hacer la union de dos grafos y asi introducimos en nuestro nuevo grafo también la union con las proteinas covid.

#protCovProtHum <- merge(x = covid, y = grafoFinal, by = c("protcov"))
union <-print_all(covid.igraph %u% rutaFinalgrafo)
plot(union)
covidprot<- covid$covid1
covidprot<- unique(covidprot)

# Podemos pintar tambien la ruta mas corta entre las del covid.Para ver la ruta directa.
grafounion <- all_shortest_paths(union, to=covidprot[1], from =covidprot[2], mode="all")

#Una vez unido formar nuevamente el dataframe dos a dos para pasarlo a grafo y pintarlo
vS3<- c()
vS4<- c()
for(i in 1:length(grafounion$res)){
  for(j in 2:length(grafounion$res[[i]])-1){
    ruta1<- grafounion$res[[i]][j]$name
    ruta2<- grafounion$res[[i]][j+1]$name
    vS3<- c(vS3, ruta1)
    vS4 <- c(vS4,ruta2 )
  }
}

grafoFinalUnion <- data.frame(vS3,vS4)
grafoFinaligraph<- graph_from_data_frame(grafoFinalUnion)
plot(grafoFinaligraph)

## average shorted past, como de cercana esta la red entre proteinas viricas.
## podemos compararlo con toda la red del proteoma para ver nuestro subgrafo como se comporta. 
## mirar componentes conexas del grafo del proteoma, coger la que es mucho mas grande que el resto.
## ver si las prot del covid caen en grafos distintos.