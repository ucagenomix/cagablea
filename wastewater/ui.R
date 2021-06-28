#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.

library(shiny)
library(leaflet)

# Define UI for application that draws a histogram
shinyUI(fillPage(
    tabsetPanel(
      tabPanel("Variants of concern", 
               pageWithSidebar(
                 titlePanel("Variants of concern data"),
                 sidebarPanel(
                   selectInput("lineage", "Lineage:", 
                               c("B.1",
                                 "P.1",
                                 "B.1.160",
                                 "B.1.351",
                                 "B.1.1.7",
                                 "A.23.1",
                                 "A.27",
                                 "B.617.2"
                               )),
                   selectInput("month", "Month:", 
                               c("oct",
                                 "nov",
                                 "dec",
                                 "jan",
                                 "feb",
                                 "mar",
                                 "apr"
                               )),
                               h3("Alterations associated with low abundant variants of concern (P.1, B.1.351, B.1.1.7, A.23.1) and choosen major lineages"),
                               p("The size of dots corresponds to total number of reads aligned to the position, the mutations on the x axis correspond to the alterations found in the samples.  The percentatage of P.1 can be estimated by the mutations that are not shared by other lineages: S:L18F, orf8:E92K show up only in February in 'La Madeleine' and 'Port'. Mutations uniquly associated with B.1.351 are visible also in February as lowfrequency alterations in Carras, Port and East Nice ")
                 ),
                 # Application title
                 mainPanel(
                   
                   plotOutput("ggplot",
                 height = "800px")
                 )
               )
             ),
      tabPanel("Lineages distribution", pageWithSidebar(
        titlePanel("lineage distribution data"),
        sidebarPanel(
          selectInput("area", "Area:", 
                      c("General",
                        "ARIANE", "BON VOYAGE", "CARABACEL", "CARRAS", "EAST GAMBETTA", "EAST JEAN MEDECIN", "EAST NICE", "HALIOTIS",
                        "LAS PLANAS", "LES MOULINS","MADELEINE","MAGNAN","MUSICIENS","NICE ETOILE","PAILLON","PORT",
                        "VIEUX NICE","WEST GAMBETTA","WEST JEAN MEDECIN","WEST NICE","FABRON","NORTH ARIANE","SOUTH ARIANE"
                      )),
                               h3("The percentage of different lineages present in waste water samples"),
                               p("The percentage was counted as the median of alteration frequencies of the reference mutations of a certain lineage found in a sample")
                 ),
        # Application title
        mainPanel(
          
          plotOutput("ggplot_2",
                 height = "800px")
        )
      )
    ),
    
    tabPanel("WW sampling stations", pageWithSidebar(
      titlePanel("WW sampling stations"),
      sidebarPanel(
        selectInput("month_ww", "Month:", 
                    c("oct", "nov", "dec", "jan", "feb", "mar", "apr")),
                    h3("The spatial distribution of the lineages present in waste water samples"),
                               p("The percentage was counted as the median of alteration frequencies of the reference mutations of a certain lineage found in a sample. Click on a barplot to get the detailed summary.")
                 ),
      # Application title
      mainPanel(
        leafletOutput("minicharts",
                 height = "800px")
      )
    )
    )
    )
  )
  
    
)
