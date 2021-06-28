#######
# data WW (ivar output) to localhost SQL DB
#05/03
# Anna
#######################################################

library(DBI)
library(RMySQL)
library(rlist)
library(stringr)
library(dplyr)

con <- DBI::dbConnect(RMySQL::MySQL(), 
                      host = "localhost",
                      user = "root",
                      dbname = "cagablea",
                      password = rstudioapi::askForPassword("Database password"))

#inserting the ivar output to created before Ww_ivar table:
#takes a path to a folder with the ivar output as an argument
#a function to create one big table from all the files in one folder:
#!!!ivar data can be added after decoding table is filled!!

ivar_data <- function(mypath, con){
  filenames <- list.files(mypath)
  #create a vector-dictionary:
  #vector with encoded proteins:
  translate_prot <- c("CDS:ENSSASP00005000003", "CDS:ENSSASP00005000002", "CDS:ENSSASP00005000004", "CDS:ENSSASP00005000006", "CDS:ENSSASP00005000010",
                      "CDS:ENSSASP00005000007", "CDS:ENSSASP00005000005","CDS:ENSSASP00005000011", "CDS:ENSSASP00005000013", "CDS:ENSSASP00005000009",
                      "CDS:ENSSASP00005000012", "CDS:ENSSASP00005000008")
  #correcponding to the codes protein names:
  names(translate_prot) <- c("ORF1ab", "ORF1ab", "S", "ORF3a", "E", "M", 'N', "ORF6", "ORF10", "ORF7a", "ORF7b", "ORF8")
  #creates a dataframe with first and last position of protein coding regions:
  prot_coord = data.frame(start = c(266,	266,	21563,	25393,	26245,	26523,	27202,	27394,	27756,	27894,	28274,	29558),
                          stop = c(13483,	21555,	25384,	26220,	26472,	27191,	27387,	27759,	27887,	28259,	29533,	29674),
                          features = c("ORF1ab",	"ORF1ab",	"S",	"ORF3a",	 "E",	 "M",	 "ORF6",	 "ORF7a",	"ORF7b",	 "ORF8",	 "N",	"ORF10"))
  i<-0
  x <- list()
  for (file in filenames){
    i<-i+1
     # (x [i]) <- paste("ivar", i, sep = "")
    print(paste(mypath, file, sep = "\\"))
    buffer_var<- read.table(paste(mypath, file, sep = "\\"), sep = '\t', header = TRUE)
    #removes the duplicated rows that ivar soetimes produces
    buffer_var <- buffer_var%>%distinct(POS, REF, ALT, REF_DP, REF_RV, REF_QUAL, ALT_DP,     
                                        ALT_RV, ALT_QUAL, ALT_FREQ, TOTAL_DP, PVAL, PASS, GFF_FEATURE, REF_CODON, 
                                        REF_AA, ALT_AA, ALT_CODON)
    buffer_var$filename <- sub(".tsv", "", file)#remove .tsv from the filename column
      #alt_dp normally corresponds to total_dp - ref_dp. In case of deletions ivar counts alt_dp of the previous base. 
    # this way the sum of alt_dp + ref_dp for deletions never corresponds to the total_dp value
    # this issue was reported on ivar github as issue #86 https://github.com/andersen-lab/ivar/issues/86
    #to change this value i count the difference between total coverage and depth of reference base:
    deletions <- which(str_detect(buffer_var$alt, "-"))
    for (i in 1:nrow(buffer_var[deletions,])){
      buffer_var$alt_dp[deletions][i] <- buffer_var$total_dp[deletions][i] - buffer_var$ref_dp[deletions][i]
      i = i + 1
    }
    #to recount alt_frequency of deletions I count the ratio of alt_dp/total_dp:
    for (i in 1:nrow(buffer_var[deletions,])){
      buffer_var$alt_freq[deletions][i] <- buffer_var$alt_dp[deletions][i]/buffer_var$total_dp[deletions][i]
      i = i + 1
    }
    #translate gff_feature which are encoded in ivar output in number format 
    
    #replace all the codes by feature from vector's names:
    buffer_var$GFF_FEATURE <- stringr::str_replace_all(buffer_var$GFF_FEATURE, setNames(names(translate_prot), translate_prot))
   
  
    #add deletions in del:28278:3 format as aa column:
    buffer_var$aa <- rep(NA, nrow(buffer_var)) 
    deletions <- which(str_detect(buffer_var$ALT, "-"))#gives a vector with the rows in samples where there are deletions
    buffer_var$aa[deletions] <- paste0("del:", buffer_var$POS[deletions], ":", nchar(buffer_var$ALT[deletions]))#pastes del: followed by the number of missing bases
    buffer_var$GFF_FEATURE[deletions] <- buffer_var$aa[deletions]
    
    #paste into aa column the features in S:T22K format, i.e. in protein coordinates:
    deletions <- which(str_detect(buffer_var$aa, "del"))#coordinates of features that excludes deletions
    i = 1
    for (i in 1:nrow(buffer_var[-deletions,])){#counts the number of codon and adds the ref amino acid and the altered amino acid
      buffer_var$aa[-deletions][i] <- (buffer_var$POS[-deletions][i] - unname(get_prot_coord[buffer_var$GFF_FEATURE[-deletions][i]]))%/%3+1
      i = i + 1
    }
    buffer_var$aa[-deletions] <- paste0("aa:", buffer_var$GFF_FEATURE[-deletions], ":", buffer_var$REF_AA[-deletions], buffer_var$aa[-deletions], buffer_var$ALT_AA[-deletions])
    
    insertStatement = DBI::sqlAppendTable(con, "ww_ivar", buffer_var, row.names = FALSE)
    x <- list.append(x, buffer_var)
    DBI::dbExecute(con, insertStatement)
    
  }
  
}
ivar_files <- list.files("D:\\data\\iVar_output\\")
#mypath <- "D:\\data\\iVar_output\\27_jan_21"
mypath <- paste0("D:\\data\\iVar_output\\", ivar_files[9])
ivar_data(mypath, con)

dbDisconnect(con)
 






