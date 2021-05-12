#################################################################
# 05.05. 2021
# Anna Diamant
# UPDATED VERSION
# builds a heatmap from VOC and variant of interest mutations 
# using pheatmap library and orders samples by month:
#################################################################
library(DBI)
library(RMySQL)
library(tidyverse)
library(rlist)
library(pheatmap)
con <- DBI::dbConnect(RMySQL::MySQL(), 
                      host = "localhost",
                      user = "root",
                      dbname = "cagablea",
                      password = "")

# the request includes the samples from all the months except the seq from 201110 and 201109 that were resequenced:
VOC <- 'SELECT * 
FROM `ww_ivar`, `pango_lineages`,`ww_barcode_decode` 
WHERE `ww_barcode_decode`.`filename` = `ww_ivar`.`filename` 
AND `ww_ivar`.`pos` = `pango_lineages`.`pos` 
AND `ww_ivar`.`alt` = `pango_lineages`.`alt`
AND `pango_lineages`.`lineage` IN ("A.27", "B.1.1", "B.1.351", "B.1.1.7", "P.1", "B.1.160",
"B.1.367", "B.1.177", "B.1.221", "B.1.474")
AND  ww_ivar.filename NOT LIKE BINARY "%201109%"
AND  ww_ivar.filename NOT LIKE BINARY "%201110%"
AND ww_barcode_decode.area NOT IN ("Lezardieux", "Pleubrian", "CergyPontoise", "Ste MAX FRAIS", "Ste MAX vieilli", "Ariane2_vieillissement", "Haliotis bis")'

# heatmap_cluster <- function(con, query){
var_jan <- dbGetQuery(con, VOC)
var_jan <- var_jan[,unique(colnames(var_jan))]#sql join is producing the dublicated columns. this lile removes them
var_jan$barcode <- sub("barcode", "", var_jan$filename)
var_jan$barcode <- sub ("cagablea_" , "" , var_jan$filename)
var_jan$barcode <- sub ("cagablea." , "" , var_jan$filename)
var_jan$barcode <-  sub("^P([A-Z 0-9]{7}_)", "", var_jan$filename)

var_jan$barcode <-  sub("210127", "jan", var_jan$barcode)
var_jan$barcode <-  sub("20210302", "feb", var_jan$barcode)
var_jan$barcode <-  sub("20201120", "nov", var_jan$barcode)
var_jan$barcode <-  sub("20201223", "dec", var_jan$barcode)
var_jan$barcode <-  sub("210224", "oct", var_jan$barcode)
var_jan$barcode <-  sub("210325", "mar", var_jan$barcode)

#var_jan$barcode <- paste(var_jan$barcode, var_jan$area, sep = "_")
#201110, 201109 are excluded in the initial request 
#lapply(var_jan$barcode, function (x) which("201110" %in% x))
#Ste MAX FRAIS, Ste MAX vieilli, 
#unique(var_jan$barcode)
var_jan$mut <- paste0(var_jan$ref,var_jan$pos, var_jan$alt)
#remove the mutations that repeat for several strains:
mut_unique <- var_jan%>%distinct(lineage, mut)%>%count(mut)%>%filter( n == 1)#mut that happen only in one lineage 
b1_mut <- var_jan%>%distinct(lineage, mut)%>%count(mut)%>%filter( n == 9)#B1 mutations 
mut_unique <- rbind(mut_unique, b1_mut)
var_jan <- var_jan[var_jan$mut%in%mut_unique$mut,]#subset the dataset to have only mutations that are unique
#create an empty data.frame in wide format:
jan_wide <-data.frame(matrix(ncol = length(unique(var_jan$mut)), nrow = length(unique(var_jan$barcode))))
colnames(jan_wide) <- unique(var_jan$mut)  #name columns as mutations to compare
row.names(jan_wide) <- unique(var_jan$barcode)
#fill up the table with the frequency data:
i = 1
while (i <= nrow(var_jan)){
  a <- var_jan[i,"alt_freq"] #read the variable in alt_freq col
  y <- var_jan$mut[i] #name of a column in jan_wide where the frequency should be written
  x <- var_jan$barcode[i]#name of a row in jan_wide where the frequency should be written
  jan_wide[c(x),c(y)] <- a #write the frequency of a mutation in a corresponding cell coded by i and b
  i = i+1
}

jan_wide[is.na(jan_wide)] <- 0 #replace NA with 0

#build a heatmap
#return(as.matrix(jan_wide))
#}

#df <-  heatmap_cluster(con, VOC)
df <- as.matrix(jan_wide)
library(pheatmap)
svg()


#my_strain_col <- var_jan[,c("mut", "lineage")]%>%distinct(mut, lineage)
#write.csv(my_strain_col, "my_strain_col.csv")#to edit it in exel
#my_strain_col%>%count(mut)%>%filter(n>1)#to decide which rows should be deleted
my_strain_col <- read.csv("my_strain_col.csv", header = T)
my_strain_col <- data.frame(my_strain_col$lineage)
rownames(my_strain_col) <- colnames(df) 
colnames(my_strain_col) <- "lineage"
rownames(my_strain_col)==my_strain_col$mut#they match

#reorder the samples' names in chronological order:
nam_ordered <- c(rownames(df)[which(stringr::str_detect(rownames(df), "oct_"))], 
                 rownames(df)[which(stringr::str_detect(rownames(df), "nov_"))],
                 rownames(df)[which(stringr::str_detect(rownames(df), "dec_"))], 
                 rownames(df)[which(stringr::str_detect(rownames(df), "jan_"))],
                 rownames(df)[which(stringr::str_detect(rownames(df), "feb_"))], 
                 rownames(df)[which(stringr::str_detect(rownames(df), "mar_"))])

df_ordered <- df[nam_ordered,]#df with samples sorted from oct to mar
rownames(df_ordered) <- nam_ordered

my_sample_col <- data.frame(month = rep(c("oct", "nov", "dec", "jan",  "feb", "mar"), c(17, 18, 19, 22, 21, 21)))
#cbind(rownames(df), my_sample_col$month)  to test if the samples names and seqeunce of months match                                        
row.names(my_sample_col) <- rownames(df_ordered)

annoCol<-list(lineage=c(A.27 = "orange", B.1.1 = "grey", B.1.1.7 = "steelblue3", B.1.160 = "indianred", B.1.177 = "palegreen3", 
                        B.1.221 = "darkblue", B.1.351 = "brown",
                        B.1.367 = "orchid4", B.1.474 = "turquoise2", P.1 = "darkolivegreen"), 
              month = c(oct = "yellow", nov = "green", dec = "blue", jan = "lightblue", feb = "orange", mar ="red"))

pheatmap(t(df_ordered), cellwidth=4, cellheight=5, fontsize= 5, cluster_cols = F, show_colnames = T, show_rownames = T,
         annotation_col = my_sample_col, 
         annotation_row =  my_strain_col, annotation_colors = annoCol)
