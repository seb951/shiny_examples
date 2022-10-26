#
library(shiny)
library(cluster)



shinyUI(fluidPage(
    
    # Application title
    titlePanel("Clustering example"),
    htmlOutput("text_cluster"),
    sidebarLayout(
        sidebarPanel(
            sliderInput("krange","krange",min = 2, max = 8, value = c(2, 5))),
            mainPanel(
                plotOutput("sil")
            )),
    
    sidebarLayout(
        sidebarPanel(
            sliderInput("k_selected","Select a K based on Silhouette plot above",min = 2, max = 8, value = 2),
            #radioButtons("pcX", "Choose two PCs:",
             #        choiceNames = paste0("PC",1:10),
              #       selected = c(1),
               #      choiceValues = 1:10
                #     )
            checkboxGroupInput("pc", "Choose two PCs",
                               choiceNames = paste0("PC",1:8),
                               selected = c(1,2),
                               choiceValues =  1:8
            ),
        ),
        mainPanel(
            plotOutput("pca")
        )),
    
    
    
    
    
)
)