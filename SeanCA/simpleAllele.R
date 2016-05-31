################################################
######Chapter 5. Migration: Spatial Models######
######Sean Wu 5/27/2016#########################
################################################


####################
###Create lattice###
####################

#A1A1=1, A2A2=2, A1A2=3
cell_init <- function(i,j,p){
  rand <- runif(1)
  if(rand < p^2){
    return(1)
  } else if(rand > 1-((1-p)^2)) {
    return(2)
  } else {
    return(3)
  }
}

cell_initV <- Vectorize(cell_init,c("i","j"))

lattice_init <- function(r,c,p){
  outer(1:r,1:c,cell_initV,p=p)
}


#########################################
###Helper functions to run simulations###
#########################################

#functions to return an index in the viable mating range, with bounding at borders to account for toroidal lattice structure
get_index <- function(i,d){
  min <- i-d
  max <- i+d
  return(floor(runif(1) * (max - min + 1)) + (min))
}

bounded_index <- function(i,d,dim){
  index <- get_index(i=i,d=d)
  if(index <= 0){
    index <- index + dim
  }
  if(index > dim){
    index <- index - dim
  }
  return(index)
}

#function to return a mate for the cell at position i,j
get_mate <- function(i,j,d,grid){
  i_prime <- bounded_index(i=i,d=d,dim=nrow(grid))
  j_prime <- bounded_index(i=j,d=d,dim=ncol(grid))
  return(grid[i_prime,j_prime])
}

#function to return offspring for two parents (A1A1=1, A2A2=2, A1A2=3)
get_offspring <- function(par1,par2){
  if(par1==1 & par2==1){ #both A1A1
    return(1)
  } else if((par1==1 & par2==3) | (par1==3 & par2==1)){ #one A1A1 one A1A2
    if(runif(1) < 0.5){
      return(1)
    } else {
      return(3)
    }
  } else if((par1==1 & par2==2) | (par1==2 & par2==1)){ #one A1A1 one A2A2
    return(3)
  } else if(par1==3 & par2==3){ #both A1A2
    rand <- runif(1)
    if(rand < 0.25){
      return(1)
    } else if(rand > 0.75){
      return(2)
    } else {
      return(3)
    }
  } else if((par1==2 & par2==3) | (par1==3 & par2==2)){ #one A1A2 one A2A2
    if(runif(1) < 0.5){
      return(2)
    } else {
      return(3)
    }
  } else { #both A2A2
    return(2)
  }
}

#function that selects a mate and produces an offspring for a single cell
cell_generation <- function(i,j,grid,d=1){
  mate <- get_mate(i=i,j=j,d=d,grid=grid)
  offspring <- get_offspring(par1=grid[i,j],par2=mate)
  return(offspring)
}

#function to run a single generation
run_generation <- function(grid){
  dim <- nrow(grid)
  iterator <- expand.grid(1:dim,1:dim) #set up iterator for mapply 
  new_grid <- mapply(cell_generation,i=iterator$Var1,j=iterator$Var2,MoreArgs=list(grid=grid,d=1),SIMPLIFY=TRUE) #use mapply to create offspring grid based on original grid
  new_grid <- matrix(data=new_grid,nrow=dim,ncol=dim) #need to turn vector into matrix; R lists matricies in column-major order
  return(new_grid)
}

run_generation <- function(grid){
  new_grid <- matrix(NA,nrow=nrow(grid),ncol=ncol(grid))
  for(i in 1:nrow(grid)){
    for(j in 1:ncol(grid)){
      new_grid[i,j] <- cell_generation(i=i,j=j,grid=grid)
    }
  }
  return(new_grid)
}


#####################
###Run Simulations###
#####################


init_grid <- lattice_init(r=100,c=100,p=0.5)
table(init_grid) #check to make sure allele frequencies make sense
sim_out <- list()

for(i in 1:100){
  if(i==1){
    grid_i <- run_generation(init_grid)
    sim_out[[i]] <- grid_i
  }
  grid_i <- run_generation(grid=grid_i)
  sim_out[[i]] <- grid_i
  print(paste("Generation i:",i))
}

image(init_grid,col=c("white","darkblue","steelblue"),axes=FALSE)
image(sim_out[[100]],col=c("white","darkblue","steelblue"),axes=FALSE)

IMAGE.dir <- "C:/Users/WuS/Desktop/migration_out/"
for(i in 1:length(sim_out)){
  write.table(x=sim_out[[i]],file=paste0(IMAGE.dir,"time",i,".csv"),col.names=FALSE,sep=",")
}
  
gg_image <- function(grid){
  melt_dat <- melt(grid)
  plot <- ggplot() +
    geom_raster(data=melt_dat,aes(x=Var1,y=Var2,fill=as.factor(value))) +
    scale_fill_manual(values=c("1"="darkblue","2"="steelblue","3"="white")) +
    guides(fill=FALSE) +
    theme(axis.text=element_blank(),axis.ticks=element_blank(),panel.background=element_blank(),panel.grid.minor=element_blank(),panel.grid.major=element_blank(),axis.title=element_blank(),title=element_blank(),panel.border=element_rect(fill=NA,colour="black",size=2))
  print(plot)
}
  
gg_image(init_grid)
gg_image(sim_out[[100]])