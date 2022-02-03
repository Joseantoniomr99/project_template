library(igraph)

grafoUnion<-read.csv("grafoUnion.csv")
Grafo<-graph_from_data_frame(grafoUnion,directed = FALSE)

#Aqui obtengo los modulos que componen el grafo, y la modularidad del grafo completo
cwt <- cluster_walktrap(Grafo)
mod <- modularity(cwt)
modularity <- modularity(Grafo,membership(cwt))

#Calculo de datos utiles en respecto a la modularidad y los modulos
communities(cwt)
sizes(cwt)
plot(wt,Grafo,vertex.label=NA)

#Posibilidad de representar dendogramas, pero la cantidad de datos hace que el resultado no sea idoneo ni util.
d <-as.dendrogram(cwt)
h <-as.hclust(cwt)