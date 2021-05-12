###########
# 8.03.2021
# Anna
# Downloading data from covidCG, nextclade etc. cleaning and adding to SQL DB
#####################################################################

library(DBI)
library(RMySQL)
library(rlist)
library(tidyverse)


#download file from covidCG:

pango <- read.csv("agg_data_lineage_All_2019-12-15-2020-12-31_top_2000.csv")
pango <- pango[,1:5]
pango_nt <- pango [,1:4]

#a function that takes as an argument the data from covidCG and returns a table to download in DB:
pango_events <- function(pango){
  pango_l <- data.frame(pos=as.character(character()),
                   ref=character(), 
                   alt=character(),
                   lineage=character(),
                   stringsAsFactors=FALSE) 
  i=0
  while (i < nrow(pango)){
    i=i+1
    if (nchar(pango[i,4]) > 0) {
      a <- pango$nt_snps [i]
      b <- sub(".", "", a)#removes first character, "["
      c <- gsub('.{1}$', '', b)#removes last character. "." stands for any, $ for end, "1" for amound of characters
      d <- strsplit(c, "[,]") #split by commas
      events <- lapply(d, function (x) strsplit(x, "[|]")) #every position with its ref and alt are stored as separate elements of a list
      j = 1
      while (j <= length(events[[1]])){
        if(events[[1]][[j]][3] == "-") {
          print(events[[1]][[j]][3])
          #delete here the first character to make it match ivar format:
          #str_sub(x,2,-1) extracts a substring from 2nd to last character:
          events[[1]][[j]][3] = paste0(events[[1]][[j]][3], str_sub(events[[1]][[j]][2],2,-1))
        }
        j <- j + 1
      }
      #append to every element of the list the lineage asociated with these mutations:
      t <- mapply(append, events[[1]], pango$lineage[i], SIMPLIFY = F)
      # print(t)
      #put every element of vector in a separate column:
      line = 1
      while (line <= length(t)){
        pango_l[nrow(pango_l)+1,] <- t[[line]]
        line = line + 1
      }
    }
    
  }
  return(pango_l)
}
temporary_var <- pango_events(pango)

con <- DBI::dbConnect(RMySQL::MySQL(), 
                      host = "localhost",
                      user = "root",
                      dbname = "cagablea",
                      password = rstudioapi::askForPassword("Database password"))

insertStatement = DBI::sqlAppendTable(con, "pango_lineages", temporary_var, row.names = FALSE)

DBI::dbExecute(con, insertStatement)

###########################################

 #download file with nextstrain clades
download.file("https://raw.githubusercontent.com/nextstrain/ncov/master/defaults/clades.tsv", "clades.tsv")
nextstrain <- read.table("clades.tsv", sep = "\t", header = T)
names(nextstrain)[3] <- "pos" #rename to make it correspond to DB structure

#define the connnection:
con <- DBI::dbConnect(RMySQL::MySQL(), 
                      host = "localhost",
                      user = "root",
                      dbname = "cagablea",
                      password = rstudioapi::askForPassword("Database password"))

insertStatement = DBI::sqlAppendTable(con, "nextstrain_clade", nextstrain, row.names = FALSE)

DBI::dbExecute(con, insertStatement)


