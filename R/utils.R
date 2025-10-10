# converts btyes to B, KB, MB or TB----
bytes.cnvrt <- function(bytes) {
  out <- ifelse(bytes < 1000, paste0(round(bytes, 2), "B"),
    ifelse(bytes < 1000000, paste0(round(bytes / 1000, 2), "KB"),
      ifelse(bytes < 1000000000, paste0(round(bytes / 1000000, 2), "MB"),
        ifelse(bytes < 1000000000000, paste0(round(bytes / 1000000000, 2), "GB"),
          paste(round(bytes / 1000000000000, 2), "TB", sep = "")
        )
      )
    )
  )
  return(out)
}


# prompt for download size----
download_data <- function(ask_before_downloading = TRUE, dl.size) {
  if (ask_before_downloading) {
    if (interactive()) {
      user_input <- readline(paste("This action will download approx ", bytes.cnvrt(dl.size), ". Would you like to proceed? (y/n): "))
      if (tolower(user_input) != "y") {
        message("Download was not performed")
        rlang::interrupt()
      }
    } else {
      message("Download was not performed, because ask_before_downloading was set to TRUE, and the environment is not interactive")
      rlang::interrupt()
    }
    message("Download beginning")
  } else {
    cat("\nTotal download approx:", bytes.cnvrt(dl.size))
    cat("\nStart time:", format(Sys.time(), "%Y-%m-%d %X"))
  }
}

# process .nc file----
nc.process <- function(nc) {
  ncfile <- ncdf4::nc_open(nc)
  nc.data <- ncdf4::ncvar_get(ncfile)
  start.day <- stringr::str_sub(string = ncfile$dim$time$units, start = nchar(ncfile$dim$time$units) - 9, end = nchar(ncfile$dim$time$units))
  dimnames(nc.data) <- list(
    ncfile$dim$lon$vals,
    ncfile$dim$lat$vals,
    as.character(as.Date(start.day) + ncfile$dim$time$vals)
  )
  ncdf4::nc_close(ncfile)
  return(nc.data)
}

# VPD function----
vpdfun <- function(tmin, tmax, relh) {
  Tave <- (tmin + tmax) / 2
  es <- 0.6108 * (exp((17.27 * Tave) / (Tave + 237.3)))
  ea <- relh / 100 * es
  VPD <- es - ea
  VPD <- VPD * 10 # Convert kpa to hpa
  return(VPD)
}

# Thermal time function----
TTfun <- function(Tci, cardT) {
  if (cardT[1] <= Tci & Tci <= cardT[2]) {
    out <- Tci
  }
  if (cardT[2] <= Tci & Tci <= cardT[3]) {
    out <- (cardT[2] / 8) * (cardT[3] - Tci)
  }
  if (Tci < cardT[1] | Tci > cardT[3]) {
    out <- 0
  }
  return(out)
}


# SMI estimation functions

# Potential Evapotranspiration----
PETthorn <- function(Tave, dl) {
  Tm <- mean(na.omit(Tave))
  N <- dl[1:length(Tave)] # day lengths
  I <- ((max(0, Tm)) / 5)^1.514 # heat index
  a <- (6.75e-07 * I^3) - (7.71e-05 * I^2) + 0.49239
  pet <- 16 * (N / 360) * ((10 * Tave) / I)^a
  pet[Tave < 0] <- 0
  return(pet)
}

# Soil moisture index function----
SMIfun <- function(prs, tmins, tmaxs, dls, AWC, cs = 2.428571, cr = 0.01236147, cd = 0.9644389, smi.start = 0.2) {
  Tave <- (tmins + tmaxs) / 2
  PETs <- PETthorn(Tave = Tave, dl = dls)
  smi <- c()
  smi[1] <- smi.start
  for (i in 1:(length(prs) - 1)) {
    smi[i + 1] <- min(1, (smi[i] / ((1 + (PETs[i] / (AWC * cs))) * cd)) + (prs[i] * cr))
  }
  return(smi)
}




#Download and extract SLGA tiffs function
dl.extrct.tifs<-function(addrs,crds){
  tmp.dir <- tempfile()
  tmp.dir <- gsub("\\", "/", tmp.dir, fixed = T)
  tmp.dir <- paste(tmp.dir, "_", 1:length(addrs), sep = "")
  options(timeout = max(50000, getOption("timeout")))
  try(utils::download.file(url = addrs, destfile = tmp.dir, method = "libcurl", quiet = T, mode = "wb"))
  finfo<-file.info(tmp.dir)
  tryagain<-which(finfo$size<1000000000 | is.na(finfo$size))
  if(length(tryagain)>0){
    try(utils::download.file(url = addrs[tryagain], destfile = tmp.dir[tryagain], method = "libcurl", quiet = T, mode = "wb"))
    }
  
  all.vals <- sapply(1:length(tmp.dir), function(d) {
    r<-terra::rast(tmp.dir[d])
    vals <- terra::extract(x = r, y = crds[,c("longitude", "latitude")], search_radius = 50)[, 2]
    trys <- 10
    rad <- 500
    while (trys > 0 & sum(is.na(vals)) > 0) { # try filling in NAS with increasing search radius
      trys <- trys - 1
      which.nas <- which(is.na(vals))
      if (length(which.nas) < 2) {
        vals[which.nas] <- unique(terra::extract(x = r, y = crds[which.nas, c("longitude", "latitude")], search_radius = rad)[, 2])
      }
      if (length(which.nas) > 1) {
        vals[which.nas] <- terra::extract(x = r, y = crds[which.nas, c("longitude", "latitude")], search_radius = rad)[, 2]
      }
      rad <- rad + 500 #increase the radius
    }
    return(vals)
  })
  file.remove(tmp.dir)
  rownames(all.vals)<-crds$Loc
  return(all.vals)
}



#extract SLGA data from API rasters function
api.extrct<-function(rasters,crds){
  AWCdata <- matrix(NA, nrow = nrow(crds), ncol = nrow(rasters), dimnames = list(
    crds$Loc,paste(paste(rasters$LowerDepth_m,rasters$UpperDepth_m,sep = "-"),"m",sep="")
  ))
  
  for (d in 1:nrow(rasters)) {
    r <- terra::rast(paste("/vsicurl/",rasters$StagingPath[d], sep = ""))
    if(nrow(crds)==1){
    vals <- unique(terra::extract(x = r, y = crds[,c("longitude", "latitude")], search_radius = 50))[, 2]
    }else{
    vals <- terra::extract(x = r, y = crds[,c("longitude", "latitude")], search_radius = 50)[, 2]
    }
    trys <- 10
    rad <- 500
    while (trys > 0 & sum(is.na(vals)) > 0) { # try filling in NAS with increasing search radius
      trys <- trys - 1
      which.nas <- which(is.na(vals))
      if (length(which.nas) < 2) {
        vals[which.nas] <- unique(terra::extract(x = r, y = crds[which.nas, c("longitude", "latitude")], search_radius = rad)[, 2])
      }
      if (length(which.nas) > 1) {
        vals[which.nas] <- terra::extract(x = r, y = crds[which.nas, c("longitude", "latitude")], search_radius = rad)[, 2]
      }
      rad <- rad + 500
    }
    AWCdata[, d] <- vals
  }
  return(AWCdata)
}


#read data from saved SLGA tiff file function
tif.read<-function(rasters,crds,tif.dir,verbose){
  files<-paste(rasters$Name,".tif",sep="")
  miss.files<-files[!files %in% dir(tif.dir)]
  if(length(miss.files)>0){ stop("SLGA TIFF files not in tif.dir:\n", paste(miss.files,collapse = "\n"),"\nTry dl.slga function again")}
  
  AWCdata<-sapply(files,function(t) {
  if (verbose == TRUE) { cat("|") }
  r <- terra::rast(paste(tif.dir,t,sep="/"))
  if(nrow(crds)==1){
    vals <- unique(terra::extract(x = r, y = crds[,c("longitude", "latitude")], search_radius = 50))[, 2]
  }else{
    vals <- terra::extract(x = r, y = crds[,c("longitude", "latitude")], search_radius = 50)[, 2]
  }
  trys <- 10
  rad <- 500
  while (trys > 0 & sum(is.na(vals)) > 0) { # try filling in NAS with increasing search radius
    trys <- trys - 1
    which.nas <- which(is.na(vals))
    if (length(which.nas) < 2) {
      vals[which.nas] <- unique(terra::extract(x = r, y = crds[which.nas, c("longitude", "latitude")], search_radius = rad)[, 2])
    }
    if (length(which.nas) > 1) {
      vals[which.nas] <- terra::extract(x = r, y = crds[which.nas, c("longitude", "latitude")], search_radius = rad)[, 2]
    }
    rad <- rad + 500
  }
  return(vals)
})
  rownames(AWCdata)<-crds$Loc
  return(AWCdata)
}





