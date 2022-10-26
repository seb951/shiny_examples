#
library(shiny)
source("~/Documents/git_repos/shiny_examples/R/phylo.R")

shinyServer(function(input, output,session) {
    output$plot1 <- renderPlot({
        plot(mtcars$wt, mtcars$mpg)
    })
    
    ###
    output$plot_phylo <- renderPlot({
        renderphylo(drop=input$dropped_tip,type= input$type)
        })
    
    
    ###
    output$coloured <-   renderPlot({
        rendercoloured(color1=input$color1,
                                      color2=input$color2,
                                      pies=input$pies,
                                      tips_colors=input$tips_colors)
    })
    
    
    ###
    output$text <- renderUI({
        para1 <- "Phylogenetic trees are one the basic tools used in a large number of evolutionary, ecological or molecular studies. In principle, they can be used to depict relationships (distance) among any kind of objects, but more usually among taxa or other biological entities (individuals, species, family, orders). In addition, having a phylogenetic tree is a necessery step in a number of analyses such as trait evolution[^1], diversification rates[^2] or molecular dating[^3]."
        para2 <- "Given their ubiquitous presence in evolutionary biology, a number of methodological approaches exist to build them (e.g. Maximum Likelihood, Neighbour Joining, Bayesian), and represent them (e.g. Figtree[^4], MEGA[^5]). Here, I will provide a few examples on how to use R to draw and manipulate them, using two common R packages (phytools[^6],also check the phytools blog[^7], and ape[^8]). I will purposefully not present any data manipulation using Hadley Wickham's ggplot2[^9] package, given that the methodology and syntax used in the tidyverse[^10] are quite different from the base R[^11] (whether it is better is another topic and a matter of opinion I think). Also note that this post was done entirely using Rmarkdown[^12] and knitr[^13], fun tools to create html, pdf, or even Word documents."
        para3 <- "I use this as a very brief introduction (and certainly non exhaustive) to plotting phylogenies in R, which I hope you will modify to fit your own datasets. In the future, I'd like to expand on this topic. Comments, questions, things you'd like to see done are welcomed!"
        HTML(paste(para1,para3, sep = '<br/><br/>'))
        
    })
    
    
})
