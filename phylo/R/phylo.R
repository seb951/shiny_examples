
#fun1
renderphylo = function(drop="",type = "phylogram"){
  #Slightly better than default plot
  if(drop!="DON'T DROP") bird.orders_dropped= drop.tip(bird.orders,drop)
  if(drop=="DON'T DROP") bird.orders_dropped= bird.orders
  if(drop=="DROP RANDOMLY 10 TIPS") bird.orders_dropped= drop.tip(bird.orders,bird.orders$tip.label[sample(1:23,10)])
  phylo2 = plot.phylo(bird.orders_dropped,type=type, edge.width = 1, label.offset = 1,cex = 0.4, main = "Figure 1 (dropping tips/type)")
  
}


#fun2
rendercoloured = function(color1="blue", color2="red",pies = "YES",tips_colors= "YES") {
    trait1 = c(sample(c(0:25),6),c(sample(25:75,6)),sample(c(75:100),11))
    names(trait1) = bird.orders$tip.label
    
    #Use contMap from phytools
    tree.trait <- contMap(bird.orders, trait1, plot = F)
    tree.trait<-setMap(tree.trait, colors=c(color1,color2))
    
    #Carefull, we now plot using the plotting tool from phytools
    plot.contMap(tree.trait, type = "phylogram", legend = 0.4*max(nodeHeights(bird.orders)), fsize = c(0.9, 0.7), lwd=3, edge.width = 2, offset = 1,mar = c(2,4,3,2))
    
    #add pies
    if(pies=="YES"){
      pies = matrix(0,nrow =bird.orders$Nnode,ncol=2)
      for(i in 1:bird.orders$Nnode) {x = sample(seq(0,1,by = 0.001),1); pies[i,] = c(x,1-x)}
      nodelabels(pie = pies, piecol = c(color1,color2), cex = 0.6)
    }
    
    if(tips_colors=="YES"){  
      tip.labels.col = c(rep(color1,11),rep(color2,12))
      tiplabels(pch=15, cex=1.2, col= tip.labels.col,offset = 0.7)
    }
    
    #Add a title manually to plot.contMap
    title("Figure 2 (adding branch colors according to a trait)",xpd = T)
  
  
}


#fun3
rendermap = function() {

renderPlot({
  
  #Matrix of random USA GPS coordinates 
  mat = matrix(c(sample(seq(from = 32, to = 45, by = 0.001),23),sample(seq(from = -120, to = -84, by = 0.001),23)),nrow = 23,ncol =2)
  rownames(mat) = bird.orders$tip.label
  
  ###World mapping function (Here we are using a USA map, but this can be changes with the database option)
  phylo.to.map(bird.orders,coords = mat, database = "usa")
  
  #Add a title manually to plot.contMap
  title("Figure 4 (geolocalised samples map)",line = -3)
},width=1400, height = 1200)

}


