library(rgdal)
library(sp)
library(maptools)
library(tlocoh)
library(foreach)
library(doParallel)
library(forecast)

#################################################

# Function to correctly format data using Kalman Smoothing to fill gaps
# Returns data frame with x and y coords in planar system, and correct POSIX datetime
# Input should be data frame with (at least) three columns with names: x, y, Datetime
format.data <- function(data, crs) {
  
  full <- data.frame(cbind(data$x, data$y, data$Datetime))
  colnames(full) <- c('long', 'lat', 'datetime')
  full$datetime <- as.POSIXct(full$datetime, tz="GMT", origin = "1970-01-01 00:00:00")
  
  long <- full$long
  z <- long
  fit <- auto.arima(long)
  kr <- KalmanSmooth(long, fit$model)
  id.na <- which(is.na(long))
  num <- ncol(kr$smooth)
  for (j in id.na) {
    z[j] <- kr$smooth[j,num]
  }
  long <- z
  
  lat <- full$lat
  z <- lat
  fit <- auto.arima(lat)
  kr <- KalmanSmooth(lat, fit$model)
  id.na <- which(is.na(lat))
  num <- ncol(kr$smooth)
  for (j in id.na) {
    z[j] <- kr$smooth[j,num]
  }
  lat <- z
  
  all <- cbind(long, lat)
  all.sp.latlong <- SpatialPoints(all[,c("long","lat")], 
                                  proj4string=CRS("+proj=longlat +ellps=WGS84"))
  all.sp.utm <- spTransform(all.sp.latlong, CRS(crs))
  all.mat.utm <- coordinates(all.sp.utm)
  all <- cbind(all.mat.utm, full$datetime)
  colnames(all) <- c('x', 'y', 'datetime')
  all <- as.data.frame(all)
  all$datetime <- as.POSIXct(all$datetime, tz="GMT", origin = "1970-01-01 00:00:00")
  
  return(all)
}

# Function for creating training-testing datasets (with approximately 80% training)
# seed argument allows for recreation of the same dataset for the sake of replication
train.test <- function (data, seed=1) {
  df <- data.frame(matrix(TRUE,nrow(data),100))
  set.seed(seed)
  
  for (i in 1:ncol(df)) {
    samp <- sample(seq(1,nrow(data),1), round(0.002222*(nrow(data))))
    for (j in 1:length(samp)) {
      for (k in 1:nrow(df)) {
        if (k == samp[j]) {
          df[k,i] <- FALSE
        }
      }
    }
  }
  
  count <- 0
  for (i in 1:ncol(df)) {
    new <- table(df[,i])[1]
    count <- count + new
  }
  
  for (i in 3:nrow(df)) {
    for (j in 1:ncol(df)) {
      if (df[i-1,j] == FALSE && df[i-2,j] != FALSE) {
        for (k in 0:98) {
          df[i+k,j] <- FALSE
        }
      }
    }
  }
  
  df <- df[1:nrow(data),]
  return(df)
}

# Function for finding test points at the center of each training-testing split
find.test.pts <- function(train.test, coords, j, crs) {
  df.temp <- data.frame(as.numeric(train.test[,j]))
  colnames(df.temp) <- "subset"
  total.pts <- sum(df.temp)
  
  # Find middle points of testing data and define as -1
  for (i in 2:nrow(df.temp)) {
    if (df.temp[i,1] == 0 && df.temp[i-1,1] != 0 && df.temp[i-1,1] != -1) {
      df.temp[i+50,1] <- -1
    }
  }
  
  df.temp <- df.temp[1:nrow(data),]
  df.temp <- cbind(coords, df.temp)
  
  # Extract the middle points for testing and record associated coordinates
  test.pts <- data.frame()
  q = 1
  for (i in 1:nrow(df.temp)) {
    if (df.temp[i,3] == -1) {
      test.pts[q,1] <- df.temp[i,1]
      test.pts[q,2] <- df.temp[i,2]
      q = q + 1
    }
  }
  
  colnames(test.pts) <- c("x", "y")
  test.pts <- SpatialPoints(test.pts[ , c("x","y")], proj4string=CRS(crs))
  
  return(test.pts)
}

# Function for performing the efficient algorithm based on:
# a) tt.split - the result of the previous function
# b) k.max - the maximum value of k that the search will reach
# c) data.lxy - a properly formatted lxy-object from tlocoh
# d) data - the original dataset (xy-coords are sufficient)
# e) crs - the a character string of the projection in proj4string form
algo.efficient <- function(tt.split, k.max, data.lxy, data, crs) {
  
  algo.grid <- function(tt.split, k.vals, s.max, num.s.vals, data.lxy, data, crs) {
    
    temp.trace <- foreach(k = 1:length(k.vals), .packages=c("tlocoh", "sp", "maptools", "forecast"), .combine='rbind') %dopar% {
      
      find.test.pts <- function(train.test, coords, j, crs) {
        df.temp <- data.frame(as.numeric(train.test[,j]))
        colnames(df.temp) <- "subset"
        total.pts <- sum(df.temp)
        
        # Find middle points of testing data and define as -1
        for (i in 2:nrow(df.temp)) {
          if (df.temp[i,1] == 0 && df.temp[i-1,1] != 0 && df.temp[i-1,1] != -1) {
            df.temp[i+50,1] <- -1
          }
        }
        
        df.temp <- df.temp[1:nrow(data),]
        df.temp <- cbind(coords, df.temp)
        
        # Extract the middle points for testing and record associated coordinates
        test.pts <- data.frame()
        q = 1
        for (i in 1:nrow(df.temp)) {
          if (df.temp[i,3] == -1) {
            test.pts[q,1] <- df.temp[i,1]
            test.pts[q,2] <- df.temp[i,2]
            q = q + 1
          }
        }
        
        colnames(test.pts) <- c("x", "y")
        test.pts <- SpatialPoints(test.pts[ , c("x","y")], proj4string=CRS(crs))
        
        return(test.pts)
      }
      
      trace <- data.frame(matrix(0,(num.s.vals+1),4))
      
      for (z in 1:(num.s.vals + 1)) {
        
        current.k_val <- k.vals[k]
        current.s_val <- (z-1) * (s.max/num.s.vals)
        
        # Calculate the nearest neighbors and create lhs object for full dataset
        full.lxy <- lxy.nn.add(data.lxy, k=current.k_val, s=current.s_val, status=F)
        full.lhs <- lxy.lhs(full.lxy, k=current.k_val, s=current.s_val, status=F)
        coords <- full.lhs[[1]]$pts@coords
        total.area <- sum(full.lhs[[1]]$hulls@data$area)
        
        # Create list for the probability values for each test/train split in df
        prob.list <- list()
        probs.log.hulls <- c()

        for (j in 1:ncol(tt.split)) {
          
          # Create a one-column data frame from df
          df1 <- tt.split[1:nrow(tt.split),j]
          
          # Create selection of hulls based on Boolean
          hulls.sel.idx <- which(df1)
          full.hulls <- hulls(full.lhs)
          selected.hulls <- full.hulls[[1]] [ full.hulls[[1]]@data$pts.idx %in% hulls.sel.idx , ]
          #temp.area <- sum(selected.hulls@data$area)
          
          # Determine the number of points in training dataset
          test.pts <- find.test.pts(tt.split, coords, j, crs)
          poly <- SpatialPolygons(selected.hulls@polygons, proj4string = CRS(crs))
          
          # Determine the number of hulls under each test point,
          overlay <- data.frame(matrix(0,length(test.pts@coords[,1]),3))
          
          for (i in 1:length(test.pts@coords[,1])) {
            overlay.list <- over(test.pts[i], poly, returnList=TRUE)
            if (length(overlay.list[[1]]) == 0) {
              overlay[i,1] <- 0 #number of hulls touching
              overlay[i,2] <- 0 #proportion of hulls being touched
              overlay[i,3] <- 1 / (total.area^2) #number of hulls over total area (prob)
            } else {
              overlay[i,1] <- length(overlay.list[[1]]) #number of hulls touching
              overlay[i,2] <- (overlay[i,1]/length(selected.hulls)) #proportion of hulls being touched
              overlay[i,3] <- overlay[i,1] / total.area #number of hulls over total area (prob)
            }
          }
          
          colnames(overlay) <- c("over", "prop", "prob.hull.area")
          log.prob.hull <- log(overlay$prob.hull.area)

          # Add values to likelihood list
          prob.list[[j]] <- as.list(overlay)
          probs.log.hulls[j] <- sum(log.prob.hull)
        }
        
        trace[z,1] <- current.s_val
        trace[z,2] <- current.k_val
        trace[z,3] <- total.area
        trace[z,4] <- sum(probs.log.hulls)
        
      } # End of s loop
      
      colnames(trace) <- c("s.val", "k.val", 'total.area',
                           'sum_log.hulls')
      
      return(trace)
      
    } # End of k loop
    
    return(temp.trace)
    
  }
  
  k.vals1 <- c()
  k.vals1[1] <- 4
  k.vals1[2:((k.max/20)+1)] <- seq(20,k.max,20)
  
  num_cores <- detectCores()
  cl<-makeCluster((num_cores - 2), outfile="")
  registerDoParallel(cl)
  
  overall.trace <- algo.grid(tt.split, k.vals = k.vals1, s.max=0.05, num.s.vals=5, 
                             data.lxy=full.lxy, data=data, crs=crs)
  
  max.id <- which.max(overall.trace$sum_log.hulls)
  max.row <- overall.trace[max.id,]
  opt.k <- max.row$k.val
  
  print(max.row)
  
  if (opt.k > 20) {
    temp.min <- opt.k - 20
  } else {
    temp.min <- 5
  }
  temp.max <- opt.k + 20
  k.vals2 <- seq(temp.min, temp.max, 5)
  
  overall.trace2 <- algo.grid(tt.split, k.vals = k.vals2, s.max=0.05, num.s.vals=5, 
                              data.lxy=full.lxy, data=data, crs=crs)
  
  max.id2 <- which.max(overall.trace2$sum_log.hulls)
  max.row2 <- overall.trace2[max.id2,]
  opt.k2 <- max.row2$k.val
  
  print(max.row2)
  
  temp.min2 <- opt.k2 - 5
  temp.max2 <- opt.k2 + 5
  k.vals3 <- seq(temp.min2, temp.max2, 2)
  
  overall.trace3 <- algo.grid(tt.split, k.vals = k.vals3, s.max=0.05, num.s.vals=50, 
                              data.lxy=full.lxy, data=data, crs=crs)
  
  stopCluster(cl)
  
  return(overall.trace3)
  
}

###############################

name.list <- c("AA", "BB", "CC", "DD")

setwd("~/")

opt.params <- data.frame(matrix(0,length(name.list),5))
for (p in 1:length(name.list)) {
  # Import data and format such that it contains x, y, and Datetime columns
  data <- read.csv(paste0(name.list[p], ".csv"))
  data$Datetime <- as.POSIXct(data$Datetime, tz="GMT", origin = "1970-01-01 00:00:00")
  crs <- '+proj=utm +south +zone=33 +ellps=WGS84'
    
  # Execute formatting command to fill gaps, then create lxy object and train/test splits
  data <- format.data(data, crs)
  full.lxy <- xyt.lxy(xy = data[,c('x','y')],
                      dt=data$datetime, 
                      id=as.character(name.list[p]), 
                      proj4string=CRS(crs),
                      dup.dt.check=FALSE)
  tt.split <- train.test(data)
  
  # Carry out efficient algorithm (three levels: 20, 5, 2)
  strt <- Sys.time()
  overall.trace <- algo.efficient(tt.split, k.max=200, 
                                  data.lxy=full.lxy, data=data, crs=crs)
  print(Sys.time() - strt)
  
  # Select optimal k and s value based on output of efficient algorithm
  opt.id <- which.max(overall.trace$sum_log.hulls)
  opt.row <- overall.trace[opt.id,]
  
  opt.params[p,1] <- name.list[p]
  opt.params[p,2:5] <- opt.row
}

