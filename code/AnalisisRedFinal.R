## Metricas de redes y nodos

library(igraph)

##Densidad de la red inicial humana:0.03
edge_density(humana.igraph, loops=F)
##Densidad de la sub red obtenida:0.002
edge_density(union,loops=F)

##Reprocidad red inicial humana:1
reciprocity(humana.igraph)
##Reprocidad de la sub red obtenida:1
reciprocity(union)
#Diametro red incial humana:5
diameter(humana.igraph, directed=F, weights = NA)
#Diametro de la sub red obtenida: 7
diameter(union, directed=F, weights = NA)
d <- get_diameter(union, directed=F)
#Representar los nodos del diametros en color gold:
#No se aprecian correctamente ya que otros nodos los cubre
vcol<- rep("gold", vcount(union))
vcol[d]<- "gray40"
ecol<- rep("gray80", ecount(union))
ecol[E(union, path=d)]<- "orange"
plot(union, vertex.color=vcol, edge.color=ecol, edge.arrow.mode=0, vertex.label=NA)

## Grado red inicial humana
degH <- degree(humana.igraph, mode="all")
max(degH) #15014
sum(degH==15014) #1
sum(degH==2) #16
min(degH) #2
hist(degH, breaks=1:vcount(humana.igraph)-1, main="Histograma del grado", xlim=c(0,15050), ylim=c(0,30))
## Grado de la sub red obtenida
deg <- degree(union, mode="all")
max(deg) #620
sum(deg==620) #1
sum(deg==1) #208
min(deg) #1
plot(union, vertex.size=deg*3, vertex.label=NA)
hist(deg, breaks=1:vcount(union)-1, main="Histograma del grado", xlim=c(0,650), ylim = c(0,20))


# Degree distribution red humana
deg.distH <- degree_distribution(humana.igraph, cumulative=T, mode="all")
plot( x=0:max(degH), y=1-deg.distH, pch=19, cex=1.2, col="orange", 
      xlab="Degree", ylab="Cumulative Frequency")
# Degree distribution red obtenida
deg.dist <- degree_distribution(union, cumulative=T, mode="all")
plot( x=0:max(deg), y=1-deg.dist, pch=19, cex=1.2, col="orange", 
      xlab="Degree", ylab="Cumulative Frequency")

#Distancia media red humana:2.040557
mean_distance(humana.igraph, directed=F)
#Distancia media red obtenida: 3.517751
mean_distance(union, directed=F)

# Asortatividad red humana: 0.1040102
assortativity_degree(humana.igraph, directed=F)
# Asortatividad red obtenida: -0.5573838
assortativity_degree(union, directed=F)




