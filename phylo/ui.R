#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
library(shiny)

library(ape)
library(phytools)
#Dataset examples
data(bird.orders)

shinyUI(fluidPage(
    
    # Application title
    titlePanel("Plotting phylogenies"),
    htmlOutput("text"),
    # Sidebar with a slider input for number of bins
        
    sidebarLayout(
        sidebarPanel(
            selectInput("dropped_tip", "Select a tip to drop",
                            choices = c("DON'T DROP","DROP RANDOMLY 10 TIPS", bird.orders$tip.label)),
            selectInput("type", "Type of phylogeny",
                        choices = c("phylogram", "cladogram", "fan", "unrooted", "radial", "tidy"))
        ),
        mainPanel(
            plotOutput("plot_phylo")
        )
    ),
    sidebarLayout(
        sidebarPanel(
            selectInput("color1", "Select a 1st color",
                        choices = rainbow(10)),
            selectInput("color2", "Select a 2nd color",
                        choices = rainbow(10)),
            selectInput("pies", "Do you want to add pie charts",
                        choices = c("YES","NO"),selected= "NO"),
            selectInput("tips_colors", "Do you want to add colors to tips",
                        choices = c("YES","NO"),selected= "NO")
        ),
        
        mainPanel(
            plotOutput("coloured")
        )
    )

)
)