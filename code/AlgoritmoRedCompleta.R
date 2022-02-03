## 0.Importación de las librerias necesarias
library(xlsx)
library(igraph)
library(dplyr)
library(gtools)
## 1.Leemos nuestros ficheros de datos:
humana <-  read.table("9606.protein.links.v11.5.txt", sep = "", header = TRUE)
covid <- read.csv("FinalNetworkTable.csv") 
## Eliminamos la columna que no necesitamos
covid$X<- NULL
##Nos quedaremos con las columnas necesarias de la tabla covid, que sera el gen covid y el gen humano en código ensamble.
humana <- humana[, 1:2]
covid <- covid[,c(2,13)]

## Tambien podemos analizar si el grafo humano es una componente conexa o esta formado por mas de una.
## Si contiene mas de una componente debemos coger la mayor, y comprobar que las humanas se encuentran en la misma componente.
components(humana.igraph)


## 2. Crearemos los objetos de tipo igraph, eliminando la direccionalidad.
humana.igraph <- graph_from_data_frame(humana)
covid.igraph <- graph_from_data_frame(covid, directed=FALSE)

## 3. Buscaremos la interaccion solo de dos proteinas unidas al covid, para optimizar tiempo y memoria
## Guardaremos los dos primeros genes humanos que interaccionan con cada proteina covid.
Bait<-c()
ensc<- c()
for (i in 1:length(covid$Bait)){
  for(z in 1: length(covid$ensc)){
    if(covid$Bait[i]== covid$Bait[z]){
      if(length(Bait[Bait==covid$Bait[z]]) <2){
        Bait<- c(Bait, covid$Bait[z])
        ensc<- c(ensc, covid$ensc[z])
      }
    }
    
  }
}
##4. Creamos el nuevo data.frame que contendra 51 proteinas humanas y todas las de covid.
covidMod<- data.frame(Bait, ensc)

##5. Buscamos todas las rutas cortas entre estas 51 proteinas humanas unidas a proteinas covid.
## Buscamos la ruta de un nodo con su siguiente.
covprot<- covidMod$ensc
rutas1<- c()
for(i in 1:length(covprot)){
    camino<- all_shortest_paths(humana.igraph, from=covprot[i], to=covprot[i+1], mode = "all") 
    rutas1 <- c(rutas1 ,camino)
}

## Nos quedamos con las rutas unicas

unicas<-c()
for (i in 1:length(rutas1)){
  u<- unique(rutas1[i]$res)
  unicas<-c(unicas,u)
}

## Crearemos dos vectores, en un bucle iremos recorriendo las rutas unicas y para cada ruta haremos una union 2 a 2.
vS1<- c()
vS2<- c()
for(i in 1:length(unicas)){
  for(j in 2:length(unicas[[i]])-1){
    ruta1<- unicas[[i]][j]$name
    ruta2<- unicas[[i]][j+1]$name
    vS1<- c(vS1, ruta1)
    vS2 <- c(vS2,ruta2 )
  }
}

## Creamos el data frame que me genera pares de dos relaciones para las proteinas humanas. 
rutaFinal<- data.frame(vS1,vS2)
rutaFinal<- unique(rutaFinal)

##Guardamos el dataframe generado para poder usarlo posteriormente sin tener que volver a compilar todo
write.csv(rutaFinal, file = "rutaFinal.csv")

rutaFinalgrafo<- graph_from_data_frame(rutaFinal, directed = FALSE)
plot(rutaFinalgrafo)

## Veremos las componentes conexas del grafo creado, y observamos que solo tenemos una componente conexa $no=1
components(rutaFinalgrafo) 

## Debemos unir los dos grafos, el que contiene las proteinas covid y ademas todas las proteinas humanas, y nuestro grafo final obtenido.
union <-print_all(covid.igraph %u% rutaFinalgrafo)
plot(union)

##Comprobemos nuevamente las componentes conexas del grafo completo, vemos que el numero de vertices aumenta porque estan 
#incluidos los faltantes, y que el número de compnentes conexas sigue siendo 1.
components(union)

## Volvemos a guardarlo con data.frame para poder usarlo más adelante.

grafoUnion <- as_long_data_frame(union)

write.csv(grafoUnion, file = "./grafoUnion.csv")

##Coloreamos el grafo dandole un color a los genes de covid
colorescluster <-  c("blue")
plot(union,
     vertex.size=degree(union)/10,
     vertex.color=colorescluster[covid$Bait],
     vertex.label= NA
     )





##OTROS MÉTODOS PROBADOS

## Usar el all shortes path para cada gen de covid con su siguiente. Tras 5 horas de compilación no obtuvimos resultados.
## La memoria ocupada era de varios gigas.
covprot<- covid$ensc
rutas<- c()
for(i in 1:length(covprot)){
    camino<- all_shortest_paths(humana.igraph, from=covprot[i], to=covprot[i+1], mode = "all") 
    rutas <- c(rutas ,camino)
  }

rutasUnicas <- unique(rutas$res)

## Usar el all shortest path para cada gen covid con todos los demas, era nuevamente inabordable.
covprot<- covid$ensc
rutas<- c()
for(i in 1:length(covprot)){
  j=i+1
  for(j in 1:length(covprot)){
    camino<- all_shortest_paths(humana.igraph, from=covprot[i], to=covprot[j], mode = "all") 
    rutas <- c(rutas ,camino)
  }
}

##Obtener combinaciones posibles de las proteinas humanas unidas al covid, 
#para posterioremente aplicar la función de buscar los caminos que las une en el proteoma (all_simple_paths()).
# Mucho tuempo de ejecucion, tras 15 horas no obenemos resultados
covprot<- covidMod$ensc
rutas<- c()
for(i in 1:length(covidMod$ensc)){
  j=i+1
  for(j in 1:length(covidMod$ensc)){
    camino<- all_simple_paths(humana.igraph, from=covidMod$ensc[i], to=covidMod$ensc[j], mode = "all", cutoff = 3) 
    rutas <- c(rutas ,camino)
  }
}

## Intento de obtener combinaciones con all_simple paths usando cuttof, nuevamente nos quedabamos sin memoria y tardaba mucho tiempo
covprot<- covidMod$ensc
rutas<- c()
for(i in 1:length(covidMod$ensc)){
    camino<- all_simple_paths(humana.igraph, from=covidMod$ensc[i], to=covidMod$ensc[i+1], mode = "all", cutoff = 3) 
    rutas <- c(rutas ,camino)
  }


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




save(rutas, file = "./rutasCreads.RData")



