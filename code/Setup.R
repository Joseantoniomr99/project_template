paquetes <- c("BiocManager","igraph","biomaRt", "clusterProfiler","CINNA","org.Hs.eg.db","xlsx","dplyr","gtools","brainGraph")
installed <- installed.packages()[,1]
if (paquetes[1] %in% installed){
  cat(paste0("El paquete ", paquetes[1]," ya está instalado\n"))
}else{
  cat(paste0("Instalando ",paquetes[1],"\n"))
  install.packages(paquetes[1])
}
for (i in 2:length(paquetes)){
  if (paquetes[i] %in% installed){
    cat(paste0("El paquete ", paquetes[i]," ya está instalado\n"))
  }else{
    cat(paste0("Instalando ",paquetes[i]))
    BiocManager::install(paquetes[i],"\n")
  }
}