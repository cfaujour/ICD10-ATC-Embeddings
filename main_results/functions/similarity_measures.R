cosine_dist <- function(x,y, sparse = T){
  if(length(which(x!=y))==0){return(1)
  }else{
    numerator <- x%*%y
    denominator <-   max(sqrt(x%*%x)*sqrt(y%*%y), 10^-9) #avoiding division by 0 for null vector  
    if(numerator/denominator < 0){result <- 0
    }else{result <- numerator/denominator}
  }
  if(sparse){result <- ifelse(result<=0.10,0, result)}
  return(result)
} 

jaccard_sim <- function(x,y){
  return(length(intersect(x,y))/(length(x)+length(y)-length(intersect(x,y))))
}

euclidean_dist <- function(x,y){
  return(sqrt(sum((x-y)^2)))
}




#Calcul de similarité à partir de classifications existantes.

#a - fonction pour obtenir l'ancetre commun le plus bas

#table : table contenant la correspondance entre code et chapitre
get_lca <- function(x,y, table){
  if(x == y){return(x)} 
  else{
    if( table[x, "CH5"] == table[y, "CH5"]){return(table[x, "CH5"])}
    else{
      if(table[x, "CH4"] == table[y, "CH4"]){return(table[x, "CH4"])}
      else{if(table[x, "CH3"] == table[y, "CH3"]){return(table[x, "CH3"])}
        else{if(table[x, "CH2"] == table[y, "CH2"]){return(table[x, "CH2"])}
          else{return("ROOT")}}
      }
    }
  }
}


# b1 fonction pour calculer la profondeur d'un terme, se base uniquement sur la nomenclature, pas besoin d'avoir les trajectoires
get_depth <- function(x,table){
  if(x == "ROOT"){return(0)}
  else{
    if(x %in% table$CH2){return(1)}
    
    else{if(x %in% table$CH3){return(2)}
      else{if(x %in% table$CH4){return(3)}
        else{if(x %in% table$CH5){return(4)}
          else{if (x %in% table$code){return(5)}
          }
        }
      }
    }
  }
}


get_depth <- function(x,table){
  
  
  col <- vector()
  for(i in 1:ncol(table)){col[i] <- x %in% table[,i]}
  col <- min(which(col))
  return(col-1)
}


# b2 fonction pour calculer l'IC d'un terme
#Nous avons besoin de la probabilité d'apparition de chaque évènement (et des sous evenements qu'il contient) dans les consommations
#Cette mesure nécessite d'avoir les trajectoires de soins


get_IC <- function(x, table, dataset){
  if (x == "ROOT"){return(0)}
  else{
    
    
    #Récupérer tous les enfants du terme x
    column <- vector()
    
    for (i in 1:ncol(table)){
      column[i] <- x %in% table[,i]
    }
    column <- min(which(column))
    sons <- table[which(table[,column] == x), column:ncol(table)]
    sons <- unique(as.vector(as.matrix(sons)))
    
    temp2 <- dataset[which(substr(dataset$var,5,1000) %in% sons),]
    
    prob <- nrow(temp2)/nrow(dataset)
    return(log(prob))}
}


#etape 2 : itérer sur l'ensemble des couples rencontrés

get_similarity_wupalmer <- function(x,y){
  temp <- get_lca(x,y , table = chapitres_CIM)
  return(2*get_depth(temp, table = chapitres_CIM)/(get_depth(x, table = chapitres_CIM) + get_depth(y, table = chapitres_CIM)))
}

get_similarity_lin <- function(x,y){
  temp <- get_lca(x,y,table = chapitres_LPP)
  return(2*get_IC(temp,subset(final_traj, traj == "ATC")) / (get_IC(x,subset(final_traj, traj == "ATC")) + get_IC(y,subset(final_traj, traj == "ATC"))))
}

#Distance de frobenius (norme 2)
dist_frobenius <- function(X,Y){A = Y - X
return(sqrt(sum(diag(t(A)%*%A)))/length(Y))}

dist_absolute <- function(X,Y){A = abs(Y-X)
return(sum(A)/length(Y))}
