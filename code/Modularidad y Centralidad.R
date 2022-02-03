library(igraph)
library(CINNA)
library(cluster)

### Centralidad
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

### Modularidad:

#Selecciono unas cuantas medidas para estudiar la centralidad:
medidasCentrality <- c("Shortest-Paths Betweenness Centrality","Page Rank","Closeness Centrality (Freeman)","Degree Centrality")
#A continuación calculo la centralidad de esas formas y ploteo un pca para ver cuál es la mejor medida de centralidad
#de las escogidas:
pdf(file = "./pcaCentralidad.pdf")
pca <- calculate_centralities(grafo,include = medidasCentrality) %>% pca_centralities(scale.unit =  TRUE, axes = c(2,2))
pca
dev.off
#La centralidad mediante cercanía de Freeman sigue una fórmula que indica el inverso de la suma de las 
#distancias, dicho de otra forma sería:
#1/(sum(i,j), i != j))
closeness <- closeness(grafo, normalized = TRUE)

#Numero de clusters con una centralidad menor de 0.1:
cat(paste0("El número de clusters con una centralidad menor que 0.1 es ", length(closeness[closeness<0.1])))

cat(paste0("El número de clusters con una centralidad esté entre 0.1 y 0.3 es ", length(closeness[closeness>=0.1 & closeness<0.3])))

cat(paste0("El número de clusters con una centralidad esté entre 0.3 y 0.5 es ",length(closeness[closeness>=0.3 & closeness<0.5])))

cat(paste0("El número de clusters con una centralidad mayor que 0.5 es ", length(closeness[closeness>=0.5])))

write.table(closeness, file = "./Closeness.csv")