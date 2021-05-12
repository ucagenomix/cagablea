#########################
# decoding file to SQL DB
#
#################

library(DBI)
library(RMySQL)

decode <- read.csv("C:\\Users\\IPMC\\Documents\\GitHub\\covid-ww-analysis\\ww_samples_analysis\\format_data_SQL\\samples_description.csv", header = T)
decode <- decode[c(1:141),c(1:3)]
names(decode) <- c("filename", "area", "month")
decode$filename[1:19] <- paste0("PAE56548_", decode$filename[1:19])

con <- DBI::dbConnect(RMySQL::MySQL(), 
                      host = "localhost",
                      user = "root",
                      dbname = "cagablea",
                      password = rstudioapi::askForPassword("Database password"))

insertStatement = DBI::sqlAppendTable(con, "ww_barcode_decode", decode, row.names = FALSE)
DBI::dbExecute(con, insertStatement)