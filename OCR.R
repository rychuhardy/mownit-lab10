library(png)
library(iterators)
library(itertools)
library(EBImage)

getPatterns <- function(filename)
{
  patterns <<- list()
  ifile <- ihasNext(ireadLines(filename))
  while(hasNext(ifile)) {
    line <- nextElem(ifile)
    
    # Convert string to list of chars 
    line <- unlist(strsplit(line, "\t"))
    
    # Second char is the represented letter
    char <- line[2]
    # From 7th character starts the matrix representation
    line <- matrix(as.numeric(line[7:length(line)]), ncol=8, nrow=16, byrow=TRUE)
    
    # Swap values (In this file black is 1 and white is 0
    # while in rasterImage it's the other way round)
    line[line==0] <- 2
    line[line==1] <- 0
    line[line==2] <- 1 
    
    # Calculate fft
    #line <- fft(line)
    
    patterns[[length(patterns)+1]] <<- line
    names(patterns)[length(patterns)] <<- char
    
  }
  
}

loadIMG <- function(filename)
{
  # Load RGB image
  img <- readPNG(filename)
  # Reduce to grayscale
  gray <- img[,,1]+img[,,2]+img[,,3]
  # Normalize
  gray <- gray/max(gray)
  # Round to 0 and 1 only (1 if gray[x] > 0.7 else 0)
  # 0.7 is chosen experimentally
  gray[gray>0.3] <- 1
  gray <- round(gray)
  # Plot result
  plot(c(0,1),c(0,1),t='n')
  rasterImage(gray, 0,0,1,1)
  return (gray)
}

plot.cluster <- function(X,C, plot=TRUE)
{
  
  
  x1 <- min(C[,1])
  y1 <- min(C[,2])
  
  x2 <- max(C[,1])
  y2 <- max(C[,2])
  
  if(plot)
    rasterImage(X[x1:x2,y1:y2], 0,0,1,1)
  
  return (X[x1:x2, y1:y2])
}

match.image <- function(image, C)
{
  #C <- dbscan(image)
  letters <- list()
  for(i in 1:length(C)) {
    

    letters[[length(letters)+1]] <- match.letter(letter)
    
  }
  return (letters)
}

match.letter <- function(letter)
{
  
  #letter <- rotate(letter, angle=270)
  
  letter <- resize(letter, w=16, h=8)
  #rasterImage(letter,0,0,1,1)
  # Calculate fft of the matched letter
  #f <- fft(letter)
  
  min <- 1
  min_val <- euc.dist(letter, patterns[[1]])
  #min_val <- Re(fft(patterns[[1]]*f,inverse=TRUE))
  #min_val <- min_val/max(min_val)
  
  for(i in 2:length(patterns)) {
    
    #val <- Re(fft(patterns[[i]]*f,inverse=TRUE))
    val <- euc.dist(letter, patterns[[i]])
    #val <- val/max(val)
    if(val > min_val) {
      min_val <- val
      min <- i
    }
  }
  return(names(patterns)[min])
    
}

# Helper function
euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))


# Adjusts clusters size based on largest cluster. This is required because the letter I
# is represented as black only matrix and thus it is impossible to match it.
resize.clusters <- function(X,C)
{ 
  # Dimenstion of each cluster squared
  diffX <- lapply(C, function(X) {return (max(X[,1])-min(X[,1]))} )
  diffY <- lapply(C, function(X) {return (max(X[,2])-min(X[,2]))} )
  
  maxX <- max(unlist(diffX))
  maxY <- max(unlist(diffY))
  # Fill smaller clusters with white points so that each has the same size
  
  for(i in 1:length(C)) {
    d <- maxX - diffX[[i]]
    
    C[[i]] <- plot.cluster(X,C[[i]],FALSE)
    
    m <- matrix(1,nrow=d/2, ncol=ncol(C[[i]]))
    C[[i]] <- rbind(m, C[[i]])
    C[[i]] <- rbind(C[[i]], m)
    
    d <- maxY - diffY[[i]]
    m <- matrix(1,ncol=d/2, nrow=nrow(C[[i]]))
    C[[i]] <- cbind(m, C[[i]])
    C[[i]] <- cbind(C[[i]], m)
  }
  return (C)
}
 

# Finds lines of text in image
find.lines <- function(X)
{
  result <- c()
  #strokewidth
  sw <- 3
  # Find top margin of a line in text
  x1 <- 1
  while(x1+sw < nrow(X)) {
  
    while(all(X[x1:(x1+sw),]==1) && (x1+sw) < nrow(X)) {
      x1 <- x1+1
    }
    x2 <- x1
    # Find bottom margin of a line in text
    while(any(X[x2:(x2+sw),]==0) && (x1+sw) < nrow(X)) {
      x2 <- x2+1
    }
    
    #Append middle of the two lines to result
    if(x1 != x2)
      result <- append(result, (x2+x1)/2)
    
    #Next line...
    x1 <- x2+1  
  }
  return (result)
}


# X is an image matrix containing zeros and ones only
#
# Eps is the maximum distance between two points of a letter.
# 1 means they are adjacent in X.
#
# Minpts is the minimum number of points to create a cluster.
#
dbscan <- function(X, eps=1, minpts=2)
{
  
  # Initialize to zeros
  Visited <- matrix(0, ncol=ncol(X), nrow=nrow(X))
  
  # Clusters
  C <- list()
    
    L <- find.lines(X)
    for(i in L) {
    #for(i in 1:nrow(X)) {
      for(j in 1:ncol(X)) {
        # Skip already visited and white points (= non letter points)
        if(Visited[i,j] == 0 && X[i,j] == 0) {
          Visited[i,j] = 1
          neighbours = regionQuery(X, i, j, eps)
          if(nrow(neighbours) >= minpts) {
            tmp <- expandCluster(X, i, j, Visited, C, neighbours, eps, minpts)
            C <- tmp$C
            Visited <- tmp$Visited
          } 
        }
      }
    }
  
  return (C)
  
}

# C is list of all current clusters
#
# neighbours is array of size n x 2 where n[k,1], n[k,2] are indexes (x,y) of a point in X
#
expandCluster <- function(X, i, j, Visited, C, neighbours, eps, minpts)
{
  # Initialize new cluster with point X[i,j]
  cluster <- matrix(c(i,j), ncol=2)
  
  # This loop doesn't work - is not updated with neighbours growing
  # And Visited array is not updated - return it from this function to dbscan
  #for(k in 1:nrow(neighbours))
  k <- 1
  while(k<= nrow(neighbours)) {
    if(Visited[neighbours[k,1], neighbours[k,2]] == 0) {
      Visited[neighbours[k,1], neighbours[k,2]] = 1
      other_neighbours <- regionQuery(X, neighbours[k,1], neighbours[k,2],eps)
      if(nrow(other_neighbours) >= minpts) {
        # join other_neighbours with neighbours (no repetitions)
        neighbours <- unique(rbind(neighbours, other_neighbours)[,1:2])
        
        # The pseudo code below isn't implemented for perfomance
        # and it doesn't matter as long as eps is 1,
        # because this would be an impossible situation.
        #
        #   if P' is not yet member of any cluster
        #     add P' to cluster C
        #
      }
    }
    k <- k+1
  }
  # Add result as a new cluster
  C[[length(C)+1]] <- neighbours
  return (list(C=C, Visited=Visited))
}


#Currrently only eps = 1 is implemented
regionQuery <- function(X, i, j, eps)
{
  result <- matrix(c(i,j), ncol=2)
  # Top left
  if(i-1 > 0 && j-1 > 0 && X[i-1, j-1]==0)
    result <- rbind(result, matrix(c(i-1,j-1), ncol=2))
  
  # Left
  if(j-1 > 0 && X[i,j-1]==0)
    result <- rbind(result, matrix(c(i,j-1), ncol=2))
  
  # Bottom left
  if(i+1 <= nrow(X) && j-1 > 0 && X[i+1,j-1]==0)
    result <- rbind(result, matrix(c(i+1,j-1), ncol=2))
  
  # Top
  if(i-1 > 0 && X[i-1,j]==0)
    result <- rbind(result, matrix(c(i-1,j),ncol=2))
  
  # Top right
  if(i-1 > 0 && j+1 <= ncol(X) && X[i-1,j+1]==0)
    result <- rbind(result, matrix(c(i-1,j+1),ncol=2))
  
  # Right
  if(j+1 <= ncol(X) && X[i,j+1]==0)
    result <- rbind(result, matrix(c(i,j+1), ncol=2))
  
  # Bottom right
  if(j+1 <= ncol(X) && i+1 <= nrow(X) && X[i+1,j+1]==0)
    result <- rbind(result, matrix(c(i+1,j+1), ncol=2))
  
  # Bottom
  if(i+1 <= nrow(X) && X[i+1,j]==0)
    result <- rbind(result, matrix(c(i+1,j),ncol=2))
  
  return(result)
}























font.size <- function(img, strokewidth=1)
{
  sw <- strokewidth
  # Find top margin of the first line of text
  x1 <- 1
  while(all(img[x1:(x1+sw),]==1)) {
    x1 <- x1+1
  }
  x2 <- x1
  # Find bottom margin of the first line of text
  while(any(img[x2:(x2+sw),]==0)) {
    x2 <- x2+1
  }
  
  # Find left margin of the first line of text
  y1 <- 1
  while(all(img[x1:x2, y1:(y1+sw)]==1)) {
    y1 <- y1+1
  }
  
  # Find right margin of first letter
  y2 <- y1
  while(any(img[x1:x2,y2:(y2+sw)]==0)) {
    y2 <- y2+1
  }
  
  rasterImage(img[x1:x2,y1:y2], 0,0,1,1)
  return (list(x1=x1, y1=y1, x2=x2, y2=y2))
  
}

find.font.size <- function(img, strokewidth = 2)
{
  # Find top margin of text
  x1 <- 1
  while(all(img[x1,]==1)) {
    x1 <- x1+1
  }
  
  #Find left margin of text
  y1 <- 1
  while(all(img[, y1]==1)) {
    y1 <- y1+1
  }
  
  # Intial position
  x <- x1
  y <- y1 
  
  # Stroke width
  sw <- strokewidth
  
  print(x)
  print(y)
  # Find the start of a letter
  while(all(img[x:(x+sw),1:y]==1) || all(img[1:x,y:(y+sw)]==1)) {
    x <- x+1
    y <- y+1
  }
  print(x)
  print(y)
  
  
  # Find the bottom of the first line of text
  while(any(img[x:(x+sw),]==0)) {
    x <- x+1
  }
  
  # Find the right side of the first letter
  while(any(img[1:x,y:(y+sw)]==0)) {
    y <- y+1
  }
  
  # Save this values for finding left and top margins of the first letter
  lx <- 0
  ly <- 0
  
  # Find top margin
  while(all(img[lx:(lx+sw), ly:y]==1)) {
    lx <- lx + 1
  }
  
  # Find left margin
  while(all(img[lx:x,ly:(ly+sw)]==1)) {
    ly <- ly + 1
  }
  
  rasterImage(img[lx:x,ly:y], 0,0,1,1)
  return (list(lx=lx, ly=ly, x=x, y=y))
}
