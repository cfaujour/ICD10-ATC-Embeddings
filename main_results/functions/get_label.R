get_label <- function(cod, ref = label_table){
  temp <- subset(ref, code == cod)
  if(nrow(temp) != 0){
    return(temp$label)}
  else{return(NA)}
}