print("Just set global variable and functions here...")

  print("Set globals")
  #load colour scales
  COLOURS_CONT <- scale_colour_gradient2()
  COLOURS_GLOBAL <-c("Blue","Yellow","Red","Green","Orange","Plum2","Purple","Brown","Cyan","Green4","Red4","Coral","Burlywood4")
  COLOURS_DIS <-scale_fill_manual(values = COLOURS_GLOBAL,aesthetics=c("fill","colour"))
  COLOURS_TWO <-scale_fill_manual(values = COLOURS_GLOBAL[c(1,5)],aesthetics=c("fill","colour"))
  MY_SEED =100
  FONTSIZE <- theme(axis.text=element_text(size=8), axis.title=element_text(size=12))
  
  
  

  
  
  
















