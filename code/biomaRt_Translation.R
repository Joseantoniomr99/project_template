library(biomaRt)
#Establecemos la conexión con la base de datos
db.Connection <- useEnsembl(biomart = "genes")
#Tomamos el dataset que queremos usar, en este caso, el de genes humanos con formato ensembl
ensembl.Human.Dataset <- useDataset(dataset = "hsapiens_gene_ensembl", mart = db.Connection)

#Cargamos el proteoma en R
interactoms <- read.table("./data/Interactoma.txt")
#Lo hemos decidido filtrar para quedarnos solo con las relaciones entre nodos cuyo score es mayor de 650
interactoma.filtered <- interactoms[interactoms$V3 > 650,]
#Quitamos el 9606. procedente de los genes humanos en String db
interactoma.filtered[,1] <- gsub("9606.",interactoma.filtered[,1],replacement = "")
interactoma.filtered[,2] <- gsub("9606.",interactoma.filtered[,2],replacement = "")

#Con la función listAttributes se ha buscado los términos claves uniprotswissprot y ensembl_peptide_id, que identifican 
#la notación de Uniprot y de Ensembl de los genes respectivamente.

#A continuación se ha realizado la query en el mart o base de datos seleccionado, de manera que queremos que nos devuelva
#de cada columna del interactoma de String db, el término de ensembl junto con el término de Uniprot correspondiente.
#Para ello, se filtra por el término de ensembl pasándole los valores de cada columna del interactoma. 

traduccion1 <- getBM(attributes=c("uniprotswissprot", "ensembl_peptide_id"), filter="ensembl_peptide_id", values=interactoma.filtered[,1], mart=ensembl.Human.Dataset)

traduccion2 <- getBM(attributes=c("uniprotswissprot", "ensembl_peptide_id"), filter="ensembl_peptide_id", values=interactoma.filtered[,2], mart=ensembl.Human.Dataset)

