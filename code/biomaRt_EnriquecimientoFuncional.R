library(clusterProfiler)
library(biomaRt)
#Establecemos la conexiÃ³n con la base de datos
db.Connection <- useEnsembl(biomart = "genes")
#Tomamos el dataset que queremos usar, en este caso, el de genes humanos con formato ensembl
ensembl.Human.Dataset <- useDataset(dataset = "hsapiens_gene_ensembl", mart = db.Connection)

#Cargamos nuestra red conexa
prots <- read.csv("./grafoUnion.csv", sep = ",")

#Nos quedamos con las columnas que tienen las proteinas
prots <- prots[,4:5]
#Recortamos el 9606. de String db que indica que son de humano
prots[,1] <- gsub(pattern = "9606.", x = prots[,1], replacement = "")
prots[,2] <- gsub(pattern = "9606.", x = prots[,2], replacement = "")
#Cogemos todas las proteinas en cada columna sin repetir
prots.Translated1 <- unique(prots[,1])
prots.Translated2 <- unique(prots[,2])

#Con la funciÃ³n listAttributes se ha buscado los tÃ©rminos claves uniprotswissprot y ensembl_peptide_id, que identifican 
#la notaciÃ³n de Uniprot y de Ensembl de los genes respectivamente.

#A continuaciÃ³n se ha realizado la query en el mart o base de datos seleccionado, de manera que queremos que nos devuelva
#de cada columna del interactoma de String db, el tÃ©rmino de ensembl junto con el tÃ©rmino de Uniprot correspondiente.
#Para ello, se filtra por el tÃ©rmino de ensembl pasÃ¡ndole los valores de cada columna del interactoma. 

#Buscamos en el mart el ENTREZID correspondiente a la notaciÃ³n ensembl de cada proteÃ­na
traduccion1 <- getBM(attributes=c("entrezgene_id", "ensembl_peptide_id"), filter="ensembl_peptide_id", values=prots.Translated1, mart=ensembl.Human.Dataset)

traduccion2 <- getBM(attributes=c("entrezgene_id", "ensembl_peptide_id"), filter="ensembl_peptide_id", values=prots.Translated2, mart=ensembl.Human.Dataset)

#Realizamos el enriquecimiento funcional de los genes de cada columna:
genes_df_GO_1 <- enrichGO(gene          = traduccion1[,1],
                        OrgDb         = org.Hs.eg.db,
                        ont           = "ALL",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.01,
                        qvalueCutoff  = 0.05,
                        readable      = TRUE)

genes_df_GO_2 <- enrichGO(gene          = traduccion2[,1],
                          OrgDb         = org.Hs.eg.db,
                          ont           = "ALL",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.01,
                          qvalueCutoff  = 0.05,
                          readable      = TRUE)


write.table(genes_df_GO_1, file = "./EnriquecimientoFuncional1.txt")

write.table(genes_df_GO_2, file = "./EnriquecimientoFuncional2.txt")

png(filename = "GenesProteinTargeting.png")
cnetplot(genes_df_GO_1)
dev.off()
png(filename = "Genes2Funciones.png")
cnetplot(genes_df_GO_2)
dev.off()
png(filename = "barplotGenes2.png")
barplot(genes_df_GO_2)
dev.off()


