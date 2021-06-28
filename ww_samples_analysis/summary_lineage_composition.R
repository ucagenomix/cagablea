####################################################
# summary on major lineages
#
####################################################


library(DBI)
library(RMySQL)
library(tidyverse)
library(rlist)
library(pheatmap)
con <- DBI::dbConnect(RMySQL::MySQL(), 
                      host = "localhost",
                      user = "root",
                      dbname = "cagablea",
                      password = rstudioapi::askForPassword("Database password"))

# the request includes the samples from all the months except the seq from 201110 and 201109 that were resequenced:
VOC <- 'SELECT * 
FROM `ww_ivar`, `pango_lineages`,`ww_barcode_decode` 
WHERE `ww_barcode_decode`.`filename` = `ww_ivar`.`filename` 
AND `ww_ivar`.`pos` = `pango_lineages`.`pos` 
AND `ww_ivar`.`alt` = `pango_lineages`.`alt` 
AND `pango_lineages`.`lineage` IN ("A.27", "B.1", "B.1.1", "B.1.351", "B.1.1.7", "P.1", "B.1.160",
"B.1.367", "B.1.177", "B.1.221", "B.1.474")
AND  ww_ivar.filename NOT LIKE BINARY "%201109%"
AND  ww_ivar.filename NOT LIKE BINARY "%201110%"
AND ww_barcode_decode.area NOT IN ("Lezardieux", "Pleubrian", "CergyPontoise", "Ste MAX FRAIS", "Ste MAX vieilli", "Ariane2_vieillissement", "Haliotis bis")
AND `ww_ivar`.`pval` < 0.00000001'
#get data:
dat <- dbGetQuery(con, VOC)
dat <- dat[,unique(colnames(dat))]#drop dublicated columns
dat$barcode <- sub("barcode", "", dat$filename)
dat$barcode <- sub ("cagablea_" , "" , dat$barcode)
dat$barcode <- sub ("cagablea." , "" , dat$barcode)
dat$barcode <-  sub("^P([A-Z 0-9]{7}_)", "", dat$barcode)#removes flowcell name
dat$barcode <- str_sub(dat$barcode, -2, -1)

dat <- dat[!duplicated(dat[,c( "filename", "pos", "ref", "alt", "ref_dp", "ref_rv", "ref_qual", "alt_dp", "alt_rv",     
                                "alt_qual", "alt_freq", "total_dp", "pval", "pass", "gff_feature", "ref_codon", "ref_aa", "alt_codon", "alt_aa",     
                               "aa", "lineage", "area", "month", "barcode")]),]
#counts the number of ref. mutations per lineage:
pango <- 'SELECT COUNT(*), `lineage` FROM `pango_lineages` 
WHERE `lineage` IN ("A.27", "B.1", "B.1.1", "B.1.351", "B.1.1.7", "P.1", "B.1.160","B.1.367", "B.1.177", "B.1.221", "B.1.474") GROUP BY `lineage`'

pango <- dbGetQuery(con, pango)

lin_stat <- dat%>%group_by( month, area, barcode,lineage)%>%summarise(median(alt_freq), min(alt_freq), max(alt_freq), sd (alt_freq))#median of alt_freq per lineage
lin_n <- dat%>%count(month, area, barcode, lineage)#number of found mut per lineage
lin_stat <- inner_join(lin_stat, lin_n, by = c("month", "area", "barcode", "lineage"))
lin_stat <- inner_join(lin_stat, pango, by = "lineage")
names(lin_stat)[10] <- "n_ref_mut"
lin_stat$percent_mut <- lin_stat$n/lin_stat$n_ref_mut
#View(lin_stat[lin_stat$percent_mut>0.8,])
lin_stat <- lin_stat[lin_stat$percent_mut>0.8,]#choose the ones there 80% or more of their ref mutations are found
lin_stat <- lin_stat%>%group_by(month, area, barcode)%>%
  mutate(lineage_percent = `median(alt_freq)`)%>%
  mutate(lineage_percent = ifelse(lineage == "B.1", NA, lineage_percent))%>%
  mutate(lineage_percent = ifelse(lineage == "B.1.1", NA, lineage_percent))%>%
  mutate(lineage_percent = ifelse(lineage == "B.1.1", `median(alt_freq)` - sum(lineage_percent, na.rm = T), lineage_percent))%>%
  mutate(lineage_percent = ifelse(lineage == "B.1",`median(alt_freq)`- sum(lineage_percent, na.rm = T), lineage_percent))

#count the minor lineages:

#Vectors with defining mutations:
#B117= c("aa:ORF1ab:T1001I", "aa:ORF1ab:A1708D","aa:ORF1ab:I2230T", "del:11288:9", "del:21765:6", "del:21991:3", "aa:S:N501Y", "aa:S:A570D", 
       # "aa:S:P681H", "aa:S:T716I", "aa:S:S982A", "aa:S:D1118H", "aa:ORF8:Q27*", "aa:ORF8:R52I", "aa:ORF8:Y73C", "aa:N:D3L", "aa:N:S235F")
b1351 = c("aa:E:P71L", "aa:N:T205I", "aa:ORF1a:K1655N", "aa:S:D80A", "aa:S:D215G", "aa:S:K417N", "aa:S:A701V", "aa:S:N501Y", "aa:S:E484K")
P1 = c("aa:ORF1ab:S1188L", "aa:ORF1ab:K1795Q", "del:11288:9", "aa:S:L18F", "aa:S:T20N", "aa:S:P26S", "aa:S:D138Y", "aa:S:R190S",
       "aa:S:K417T", "aa:S:E484K", "aa:S:N501Y", "aa:S:H655Y", "aa:S:T1027I", "aa:ORF3a:G174C", "aa:ORF8:E92K", 
       "aa:N:P80R")

A.23.1 <- c("aa:S:F157L", "aa:S:V367F", "aa:S:Q613H", "aa:S:P681R")

B.1.525 <- c("aa:orf1ab:L4715F", "aa:S:Q52R", "aa:S:E484K", "aa:S:Q677H", "aa:S:F888L", "aa:E:L21F", "aa:E:I82T", "del:11288:9", "del:21765:6", "del:28278:3")
#Looking for >50% of the 9 defining B.1.351 SNPs - defined in b1351 vector:
dat_b1351 <- distinct(dat, aa, .keep_all = T)%>%group_by( month, area, barcode)%>%filter(aa%in%b1351)%>%
  summarise(n(), median(alt_freq), min(alt_freq), max(alt_freq), sd (alt_freq))%>%filter(`n()`>4)
#Looking for at least 50% of the 16 defining P1 SNPs - defined in P1 vector 
dat_p1 <- distinct(dat, month, area, barcode, pos, ref, alt, aa, .keep_all = T)%>%group_by( month, area, barcode)%>%filter(aa%in%P1)%>%
  summarise(n(), median(alt_freq), min(alt_freq), max(alt_freq), sd (alt_freq))%>%filter(`n()`>1)
#Looking for any defining A.23.1 SNPs - defined in A.23.1 vector 
dat%>%group_by(month, area, barcode)%>%filter(aa%in%A.23.1)%>%
  summarise(n(), median(alt_freq), min(alt_freq), max(alt_freq), sd (alt_freq))%>%filter(`n()`>1)#nothing
#Looking for defining B.1.525 SNPs - defined in B.1.525 vector
dat%>%group_by(month, area, barcode)%>%filter(aa%in%B.1.525)%>%
  summarise(n(), median(alt_freq), min(alt_freq), max(alt_freq), sd (alt_freq))%>%filter(`n()`>4)#nothing 

#add summaries on minor lineages to lin_stat table:


write.csv( rbind(lin_stat2, lin_stat3),"lin_stat1.csv")



