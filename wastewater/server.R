#Anna Diamant
#May 2021
#This file defines the functions needed to vizualize the function that ui.r is calling

library(shiny)
library(ggmap)
library(ggplot2)
library(rgdal)
library(leaflet)
library(leaflet.extras)
library(rlist)
library(purrr)
library(dplyr)
library(sp)
library(raster)
library(leaflet)
library(rgdal)
library(rgeos)
library(stats)
library(ggplot2)
library(leafpop)
library(leaflet.minicharts)
library(RMySQL)
library(stringr)
# con is a variable needed to establish connection to the database. Add the password and the
# username to run the app:
# con <- DBI::dbConnect(RMySQL::MySQL(),
#                       host = "caire.ipmc.cnrs.fr",
#                       user = "....",
#                       dbname = "cagablea",
#                       port = 3306,
#                       password = "....")
area_rendered <- FALSE
plot_lineage_per <- function(area){
    if (area_rendered == FALSE) {
    VOC <- 'SELECT pango_lineages.lineage, ww_ivar.alt_freq, ww_ivar.pos, ww_ivar.alt, ww_ivar.ref, ww_ivar.aa, ww_barcode_decode.month, ww_barcode_decode.area, ww_barcode_decode.filename AS filename
FROM `ww_ivar`, `pango_lineages`,`ww_barcode_decode` 
WHERE `ww_barcode_decode`.`id` = `ww_ivar`.`barcode_id` 
AND `ww_ivar`.`pos` = `pango_lineages`.`pos` 
AND `ww_ivar`.`alt` = `pango_lineages`.`alt`
AND `pango_lineages`.`lineage` IN ("A.27", "B.1", "B.1.1", "B.1.351", "B.1.1.7", "P.1", "B.1.160",
"B.1.367", "B.1.177", "B.1.221", "B.1.474")
AND ww_barcode_decode.month IN ("oct", "nov", "dec", "jan", "feb", "mar", "apr")'
    #get data:
    dat <- dbGetQuery(con, VOC)
    dat <- dat[,unique(colnames(dat))]#drop dublicated columns
    dat$barcode <- str_sub(dat$filename, -2, -1)
    
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
    A.23.1 <- c("aa:S:F157L", "aa:S:V367F", "aa:S:Q613H", "aa:S:P681R")
    b1351 = c("aa:E:P71L", "aa:N:T205I", "aa:ORF1a:K1655N", "aa:S:D80A", "aa:S:D215G", "aa:S:K417N", "aa:S:A701V", "aa:S:N501Y", "aa:S:E484K")
    P1 = c("aa:ORF1ab:S1188L", "aa:ORF1ab:K1795Q", "del:11288:9", "aa:S:L18F", "aa:S:T20N", "aa:S:P26S", "aa:S:D138Y", "aa:S:R190S", "aa:S:K417T", "aa:S:E484K", "aa:S:N501Y", "aa:S:H655Y", "aa:S:T1027I", "aa:ORF3a:G174C", "aa:ORF8:E92K", "aa:N:P80R")
    B.1.525 <- c("aa:orf1ab:L4715F", "aa:S:Q52R", "aa:S:E484K", "aa:S:Q677H", "aa:S:F888L", "aa:E:L21F", "aa:E:I82T", "del:11288:9", "del:21765:6", "del:28278:3")

      
    #Pangolin assigns B.1.351 to sequences with at least 5 of the 9 defining B.1.351 SNPs - defined in b1351 vector:
    dat%>%group_by( month, area, barcode)%>%filter(aa%in%b1351)%>%
        summarise(n(), median(alt_freq), min(alt_freq), max(alt_freq), sd (alt_freq))%>%filter(`n()`>4)
    #Pangolin assigns P1 to sequences with at least 1o of the 9 defining P1 SNPs - defined in p1 vector 
    dat%>%group_by( month, area, barcode)%>%filter(aa%in%P1)%>%
        summarise(n(), median(alt_freq), min(alt_freq), max(alt_freq), sd (alt_freq))%>%filter(`n()`>7)
    #Pangolin assigns A.23.1 to sequences with at least defining A.23.1 SNPs - defined in A.23.1 vector 
    dat%>%group_by(month, area, barcode)%>%filter(aa%in%A.23.1)%>%
        summarise(n(), median(alt_freq), min(alt_freq), max(alt_freq), sd (alt_freq))%>%filter(`n()`>1)#nothing
    #Pangolin assigns B.1.525 to sequences with at least defining B.1.525 SNPs - defined in B.1.525 vector
    dat%>%group_by(month, area, barcode)%>%filter(aa%in%B.1.525)%>%
        summarise(n(), median(alt_freq), min(alt_freq), max(alt_freq), sd (alt_freq))%>%filter(`n()`>1)#nothing 
    #correct B.1. percentage:
    lin_stat[101, c("lineage_percent")] <- 0.1
    lin_stat[which(lin_stat$lineage_percent<0),]$lineage_percent <- 0
    
    lin_stat$month <- factor(lin_stat$month, levels = c("oct", "nov", "dec", "jan", "feb", "mar", "apr"))
    a <- lin_stat[lin_stat$area== area, ]
    
    area_rendered <- TRUE
    }
    if (area == 'General') {
        return (ggplot(lin_stat, aes(month, `median(alt_freq)`, fill = lineage))+
                geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
                #geom_text(aes(label=lineage), position=position_dodge(width=0.9), vjust=-0.25, size = 3.5)+
                facet_wrap(~area)+
                theme(axis.text.x = element_text(angle = 90, size = 14)))
    }
    return (ggplot(a, aes(month, lineage_percent, fill = lineage))+
        geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
        facet_wrap(~area)+
        theme(axis.text.x = element_text(angle = 90, size = 14)))
}
#This part defines the view of the page with bubble plots:
bubbless <-  function (lineage, mon){
        dat <- 'SELECT * FROM `ww_ivar`, `ww_barcode_decode`, `pango_lineages` WHERE `ww_ivar`.`barcode_id`= `ww_barcode_decode`.`id` AND `ww_barcode_decode`.`month` IN ("oct", "nov", "dec", "jan", "feb", "mar", "apr") AND `ww_ivar`.`pos`=`pango_lineages`.`pos` AND `ww_ivar`.`alt`= `pango_lineages`.`alt` AND `pango_lineages`.`lineage` = \''
        dat <- paste0(dat, lineage)
        dat <- paste0(dat, '\' AND ww_barcode_decode.month = \'')
        dat <- paste0(dat, mon)
        dat <- paste0(dat, '\'')
        dat <-  dbGetQuery(con, dat)
        dat$mut <-  paste(dat$ref, dat$pos, dat$alt)
        dat$aa <- gsub("aa:", "", dat$aa)

        dat$barcode <- paste(dat$month, dat$area, sep = "_")
        return (ggplot(dat[,c("month", "area","barcode", "alt_freq", "total_dp", "mut", "aa")], aes(mut, alt_freq, size = total_dp))+
                    geom_point(alpha=0.7, colour = "yellow")+
                    theme(axis.text.x = element_text(angle = 90))+
                    geom_text(aes(label = aa), size = 4, check_overlap = T, nudge_x = 0.05, colour = "black")+
                    theme(axis.text=element_text(size=14), axis.title=element_text(size=14,face="bold"))+
                    labs(x = "alterations", y = "alteration frequencies")+
                    facet_wrap(~area))
}
#this part defines the view of the page with a map and minicharts on it (barplots):
map_with_charts <- function(mon) {
    VOC <- 'SELECT pango_lineages.lineage, ww_ivar.alt_freq, ww_ivar.pos, ww_ivar.alt, ww_ivar.ref, ww_ivar.aa, ww_barcode_decode.month, ww_barcode_decode.area, ww_barcode_decode.filename AS filename
FROM `ww_ivar`, `pango_lineages`,`ww_barcode_decode` 
WHERE `ww_barcode_decode`.`id` = `ww_ivar`.`barcode_id` 
AND `ww_ivar`.`pos` = `pango_lineages`.`pos` 
AND `ww_ivar`.`alt` = `pango_lineages`.`alt`
AND `pango_lineages`.`lineage` IN ("A.27", "B.1", "B.1.1", "B.1.351", "B.1.1.7", "P.1", "B.1.160",
"B.1.367", "B.1.177", "B.1.221", "B.1.474")
AND ww_barcode_decode.month IN ("oct", "nov", "dec", "jan", "feb", "mar", "apr")'
    #get data:
    dat <- dbGetQuery(con, VOC)
    dat <- dat[,unique(colnames(dat))]#drop dublicated columns
    dat$barcode <- str_sub(dat$filename, -2, -1)
    
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
    lin_stat <- lin_stat[lin_stat$percent_mut>0.7,]#choose the ones there 80% or more of their ref mutations are found
    lin_stat <- lin_stat%>%group_by(month, area, barcode)%>%
        mutate(lineage_percent = `median(alt_freq)`)%>%
        mutate(lineage_percent = ifelse(lineage == "B.1", NA, lineage_percent))%>%
        mutate(lineage_percent = ifelse(lineage == "B.1.1", NA, lineage_percent))%>%
        mutate(lineage_percent = ifelse(lineage == "B.1.1", `median(alt_freq)` - sum(lineage_percent, na.rm = T), lineage_percent))%>%
        mutate(lineage_percent = ifelse(lineage == "B.1",`median(alt_freq)`- sum(lineage_percent, na.rm = T), lineage_percent))
    #correct B.1. percentage:
    lin_stat[101, c("lineage_percent")] <- 0.1 #recount the percentage for B.1.1 in presence of A.27
    lin_stat[which(lin_stat$lineage_percent<0),]$lineage_percent <- 0 #to remove negative percentages that appear by absence of lineage in question
    
    ################################################################################################################################################
    #turn lin_stat table into a matrix that can be plotted with minicharts:
    make_wide_df <- function(long_data){
        #create an empty data.frame in wide format:
        wide <-data.frame(matrix(ncol = length(unique(long_data$lineage)), nrow = length(unique(long_data$area))))
        colnames(wide) <- unique(long_data$lineage)  #name columns as mutations to compare
        row.names(wide) <- unique(long_data$area)
        #fill up the table with the frequency data:
        i = 1
        while (i <= nrow(long_data)){
            a <- long_data[i,"lineage_percent"] #read the variable in alt_freq col
            y <- long_data$lineage[i] #name of a column in wide where the frequency should be written
            x <- long_data$area[i]#name of a row in wide where the percentage should be written
            wide[c(x),c(y)] <- a #write the percentage in corresponding cell coded by i and b
            i = i+1
        }
        
        wide[is.na(wide)] <- 0 #replace NA with 0
        return(as.matrix(wide)) 
    }
    
    x <- make_wide_df(lin_stat[lin_stat$month==mon,])
    x <- x[ order(row.names(x)), ]
    markers <- read.csv("/data/data_diamant/minicharts_map/markers_minicharts.csv", header = T)
    markers <- markers[1:23, 1:4]
    markers <- markers[order(markers$name),]
    return (leaflet()%>%
        addTiles()%>%
        addLabelOnlyMarkers(lng = markers[markers$name%in%rownames(x),]$lon, 
                            lat = markers[markers$name%in%rownames(x),]$lat, 
                            label = markers[markers$name%in%rownames(x),]$name, 
                            labelOptions = labelOptions(noHide = T, direction = "side", textOnly = T, textsize = "14px"))%>%
        addMinicharts(
            markers[markers$name%in%rownames(x),]$lon, markers[markers$name%in%rownames(x),]$lat,
            chartdata = x,
            type = "bar",
            opacity = 1,
            width = 100,
            height = 100,
            legend = T,
            showLabels = T
        ))
    
}
# Define server logic required to draw a histogram
shinyServer(function(input, output) {
    lineage <- reactive({input$lineage})
    month <- reactive({input$month})
    month_ww <- reactive({input$month_ww})
    area <- reactive({input$area})
    output$ggplot <- renderPlot({
        bubbless(lineage(), month())
    })
    output$ggplot_2 <- renderPlot({
        plot_lineage_per(area())
    })
    output$minicharts <- renderLeaflet({
        map_with_charts(month_ww())
    })

})

