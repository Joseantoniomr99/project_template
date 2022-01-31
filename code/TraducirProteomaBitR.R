#### Filtrado de 950+ score

redProteoma <- read.table("./HumanProteome.txt", sep = "", header = TRUE)

redProt.Filtered <- redProteoma[redProteoma[,3]>950,]

redProt.Filtered[,1] <- gsub("9606.", redProt.Filtered[,1], replacement = "")

redProt.Filtered[,2] <- gsub("9606.", redProt.Filtered[,2], replacement = "")

red1.Traduccion <- bitr(redProt.Filtered[,1],fromType = "ENSEMBLPROT", toType = "SYMBOL", OrgDb = "org.Hs.eg.db")

red2.Traduccion <- bitr(redProt.Filtered[,2],fromType = "ENSEMBLPROT", toType = "SYMBOL", OrgDb = "org.Hs.eg.db")

cnt <- 1
for (i in 1:length(redProt.Filtered[,1])){
  oldNode1 <- redProt.Filtered[i,1]
  newNode1 <- red1.Traduccion[red1.Traduccion[,1] == oldNode1,][,2]
  oldNode2 <- redProt.Filtered[i,2]
  newNode2 <- red2.Traduccion[red2.Traduccion[,1] == oldNode2,][,2]
  if(length(newNode1>0)) {
    redProt.Filtered[i,1] <- newNode1[1]
  }
  if (length(newNode2>0)){
    redProt.Filtered[i,2] <- newNode2[1]
  }
  cnt <- cnt + 1
  
  
}
i <- 1
while (i < length(redProt.Filtered[,1])){
  if (substr(redProt.Filtered[i,1],1,4) == "ENSP"){
    redProt.Filtered <- redProt.Filtered[-i,]
  }
  else if (substr(redProt.Filtered[i,2],1,4) == "ENSP"){
    redProt.Filtered <- redProt.Filtered[-i,]
  }
  else{
    i <- i+1
  }
}

red950 <- redProt.Filtered

write.csv(redProt.Filtered, file = "./ProteomaTraducido.csv")

save(redProt.Filtered, file = "./ProteomaTraducido.RData")


#### Filtrado de 800 a 950 score

redProteoma <- read.table("9606.protein.links.v11.5.txt", sep = "", header = TRUE)

redProt.Filtered <- redProteoma[redProteoma[,3]<=950&redProteoma[,3]>=800,]

redProt.Filtered[,1] <- gsub("9606.", redProt.Filtered[,1], replacement = "")

redProt.Filtered[,2] <- gsub("9606.", redProt.Filtered[,2], replacement = "")

red1.Traduccion <- bitr(redProt.Filtered[,1],fromType = "ENSEMBLPROT", toType = "SYMBOL", OrgDb = "org.Hs.eg.db")

red2.Traduccion <- bitr(redProt.Filtered[,2],fromType = "ENSEMBLPROT", toType = "SYMBOL", OrgDb = "org.Hs.eg.db")

cnt <- 1
for (i in 1:length(redProt.Filtered[,1])){
  oldNode1 <- redProt.Filtered[i,1]
  newNode1 <- red1.Traduccion[red1.Traduccion[,1] == oldNode1,][,2]
  oldNode2 <- redProt.Filtered[i,2]
  newNode2 <- red2.Traduccion[red2.Traduccion[,1] == oldNode2,][,2]
  if(length(newNode1>0)) {
    redProt.Filtered[i,1] <- newNode1[1]
  }
  if (length(newNode2>0)){
    redProt.Filtered[i,2] <- newNode2[1]
  }
  cnt <- cnt + 1
  
  
}
i <- 1
while (i < length(redProt.Filtered[,1])){
  if (substr(redProt.Filtered[i,1],1,4) == "ENSP"){
    redProt.Filtered <- redProt.Filtered[-i,]
  }
  else if (substr(redProt.Filtered[i,2],1,4) == "ENSP"){
    redProt.Filtered <- redProt.Filtered[-i,]
  }
  else{
    i <- i+1
  }
}

red800 <- redProt.Filtered

write.csv(redProt.Filtered, file = "./ProteomaTraducido(800-950).csv")

save(redProt.Filtered, file = "./ProteomaTraducido(800-950).RData")


#### Filtrado de 650 a 800 score

redProteoma <- read.table("HumanProteinLinks.txt", sep = "", header = TRUE)

redProt.Filtered <- redProteoma[redProteoma[,3]<=800&redProteoma[,3]>=650,]

redProt.Filtered[,1] <- gsub("9606.", redProt.Filtered[,1], replacement = "")

redProt.Filtered[,2] <- gsub("9606.", redProt.Filtered[,2], replacement = "")

red1.Traduccion <- bitr(redProt.Filtered[,1],fromType = "ENSEMBLPROT", toType = "SYMBOL", OrgDb = "org.Hs.eg.db")

red2.Traduccion <- bitr(redProt.Filtered[,2],fromType = "ENSEMBLPROT", toType = "SYMBOL", OrgDb = "org.Hs.eg.db")

cnt <- 1
for (i in 1:length(redProt.Filtered[,1])){
  oldNode1 <- redProt.Filtered[i,1]
  newNode1 <- red1.Traduccion[red1.Traduccion[,1] == oldNode1,][,2]
  oldNode2 <- redProt.Filtered[i,2]
  newNode2 <- red2.Traduccion[red2.Traduccion[,1] == oldNode2,][,2]
  if(length(newNode1>0)) {
    redProt.Filtered[i,1] <- newNode1[1]
  }
  if (length(newNode2>0)){
    redProt.Filtered[i,2] <- newNode2[1]
  }
  cnt <- cnt + 1
  
  
}
i <- 1
while (i < length(redProt.Filtered[,1])){
  if (substr(redProt.Filtered[i,1],1,4) == "ENSP"){
    redProt.Filtered <- redProt.Filtered[-i,]
  }
  else if (substr(redProt.Filtered[i,2],1,4) == "ENSP"){
    redProt.Filtered <- redProt.Filtered[-i,]
  }
  else{
    i <- i+1
  }
}

red650 <- redProt.Filtered

write.csv(redProt.Filtered, file = "./ProteomaTraducido(650-800).csv")

save(redProt.Filtered, file = "./ProteomaTraducido(650-800).RData")


#### Union Redes 650+

TempData1 <- rbind(red650, red800)

ProteomaTCompleto <- rbind(TempData1, red950)


write.csv(ProteomaTCompleto, file = "./ProteomaBitR650+.csv")

save(ProteomaTCompleto, file = "./ProteomaBitR650+.RData")
