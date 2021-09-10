PCbiplot <- function(PC,x = 'PC1', y = 'PC2',values = TRUE) {
  # PC being a prcomp object
  
  x1 <- rlang::sym(x)
  y1 <- rlang::sym(y)
  PC_mod <- PC$model[[1]]
  data <- data.frame(obsnames=PC$qualifiers[[1]]$pl_code, PC_mod$x)
  plot <- ggplot(data, aes(x = !!x1, y = !!y1)) + 
    geom_hline(aes(yintercept = 0), size=.2) + 
    geom_vline(aes(xintercept = 0), size=.2)
  
  if(values){
    plot <- plot + geom_text(alpha=.2, size=3, aes(label=obsnames)) 
  }
  
  datapc <- data.frame(varnames=rownames(PC$model[[1]]$rotation), PC$model[[1]]$rotation)
  mult <- min(
    (max(data[,y]) - min(data[,y])/(max(datapc[,y])-min(datapc[,y]))),
    (max(data[,x]) - min(data[,x])/(max(datapc[,x])-min(datapc[,x])))
  )
  datapc <- transform(datapc,
                      v1 = 1.2 * mult * (get(x)),
                      v2 = 1.2 * mult * (get(y)),
                      v1_1 = 1 * mult * (get(x)),
                      v2_1 = 1 * mult * (get(y)) 
  )
  plot <- plot + coord_equal() + 
    geom_text(data=datapc, aes(x=v1, y=v2, label=varnames), size = 5, vjust=1, color="red")
  plot <- plot + geom_segment(data=datapc, aes(x=0, y=0, xend=v1_1, yend=v2_1), arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="red")
  plot
}