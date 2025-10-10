#' @title Get weather data from BARRA-R2
#'
#' @description Extract weather data for Australia from the [BARRA-R2](https://opus.nci.org.au/spaces/NDP/pages/264241166/BOM+BARRA2+ob53)
#' weather data resource for a set of environments with defined latitude and longitude coordinates. BARRA-R2 data runs from Jan 1979 to Sept 2024.
#'
#' @param Envs Vector of environment names character strings.
#' @param Lats Vector of latitude numeric values for each environment in the same order as `Envs`.
#' @param Lons Vector of longitude numeric values for each environment in the same order as `Envs`.
#' @param plus.yr Logical. Should the subsequent years weather data also be downloaded? This may be needed if the estimated crop growth stages in the `get.W.ECs()`
#' function extend after the same year of sowing. 
#' @param Years Vector of year integer values for each environment in the same order as `Envs`.
#' @param ncores Number (integer) of cores to use for parallel processing of gridded data over muliple years. Use `1` to run in series. The default (`NULL`) will
#' use the maximum available cores up to 5. If running in parallel, an output log text file will be created in the working directory.
#' @param verbose Logical. Should progress be printed? Default = TRUE.
#' @param dlprompt Logical. Should the user be prompted approve the total download size? Default = TRUE.
#'
#' @details
#' Weather variables returned include:
#' * `daily_rain` - Daily rainfall (mm)
#' * `max_temp` - Maximum temperature (°C)
#' * `min_temp` - Minimum temperature (°C)
#' * `vp_deficit` - Vapour pressure deficit (hPa)
#' * `radiation` - Solar exposure, consisting of both direct and diffuse components (MJ m<sup>-2</sup>)
#' * `day_lengths` - Time between sunrise and sunset (h) not taken from BARRA-R2
#' 
#' VPD in hPa is calculated as \eqn{ VPD = 10(es - ea) }, where \deqn{ es = 0.6108 \times \exp(\frac{17.27 \times T_{ave}}{T_{ave} + 237.3}) },
#' \eqn{T_{ave} } is the mean temperature in °C, \deqn{ ea = \frac{RH}{100} \times es }, and \eqn{RH} is the relative humidity (%).
#' 
#'
#' @returns A list of length 2:
#' * `$data` is a list of matrices of weather data for each weather variable.
#' Each data matrix has environment names as rows and days of the year as columns
#' * `$Env.info` is a data frame of environment names and coordinate values for environments included in the data.
#'
#' @seealso [get.SILO.weather()]
#'
#' @references
#' Su, C.H., Dharssi, I., Le Marshall, J., Le, T., Rennie, S., Smith, A., Stassen, C., Steinle, P., Torrance, J., Wang, C. and Warren, R.A., 2022.
#'     [BARRA2: Development of the next-generation Australian regional atmospheric reanalysis](http://www.bom.gov.au/research/publications/researchreports/BRR-067.pdf).
#'      \emph{Bureau of Meteorology}.
#'
#'
#' @export

get.BARRA.weather <- function(Envs,
                              Lats,
                              Lons,
                              Years,
                              plus.yr = FALSE,
                              ncores = NULL,
                              verbose = TRUE,
                              dlprompt = FALSE) {
  Years <- as.integer(as.character(Years))
  years <- unique(Years)
  if(plus.yr) years <- unique(c(years,years+1))
  Envs <- as.character(Envs)
  var.units <- c()
  mons <- stringr::str_pad(1:12, 2, pad = "0")
  vars <- c("pr", "tasmax", "tasmin", "hurs", "rsdt")

  # Check for errors
  if (length(unique(c(length(Envs), length(Lats), length(Lons), length(Years)))) > 1) {
    print(sapply(list("Envs" = Envs, "Lats" = Lats, "Lons" = Lons, "Years" = Years), length))
    stop("Lengths of Envs, Lats, Lons or Years differ")
  }
  if (!is.numeric(Lats)) stop("Lat values not numeric")
  if (!is.numeric(Lons)) stop("Lon values not numeric")
  if (sum(duplicated(Envs)) > 0) stop(paste("Duplicated Envs:", Envs[duplicated(Envs)]))
  if(plus.yr & sum(!years %in% 1979:2024) > 0) cat(crayon::yellow("Additional year of weather data outside of range (1979:2024).\n"))
  if (sum(!years %in% 1979:2023) > 0) stop("Years out of range of BARRA R2 data (Jan 1979 to Sept 2024)")
  if (sum(Lons < 88.48 | Lons > 207.39) > 0) stop("Lon out of range of BARRA data: 88.48 to 207.39")
  if (sum(Lats < -57.97 | Lats > 12.98) > 0) stop("Lats out of range of BARRA data: -57.97 to -12.98")
  if(!capabilities("libcurl")){warning("libcurl is not supported!")}

  dl.size <- 40000000 * length(vars) * length(Years) * length(mons)
  if (verbose) download_data(dlprompt, dl.size)

  if (is.null(ncores)) {
    ncores <- min(parallel::detectCores(), length(Years))
  }
  ncores <- min(ncores, length(Years))

  if (verbose) { cat("\nDownloading .nc files...\n") }
  #DL files in series
  all.vars.weather <- list()
  for (v in seq_along(vars)) {
    if (verbose) cat(vars[v], ":", sep = "")
  
    for (y in seq_along(years)) {
      if (verbose) cat(years[y], "|", sep = "")
      
      addrs <- paste("https://thredds.nci.org.au/thredds/fileServer/ob53/output/reanalysis/AUS-11/BOM/ERA5/historical/hres/BARRA-R2/v1/day/",
                 vars[v], "/latest/", vars[v], "_AUS-11_ERA5_historical_hres_BOM_BARRA-R2_v1_day_", years, mons, "-", years, mons, ".nc",
                 sep = ""
                 )
      tmp.dir <- tempdir()
      tmp.dir <- gsub("\\", "/", tmp.dir, fixed = T)
      tmp.files <- paste(tmp.dir,"/",paste(vars[v], years[y], mons, sep="_"), ".nc", sep = "")
      options(timeout = max(80000, getOption("timeout")))
      try(utils::download.file(url = addrs, destfile = tmp.files, method = "libcurl", quiet = T, mode = "wb"))
      
      finfo<-file.info(tmp.files)
      tryagain<-which(finfo$size<30000000 | is.na(finfo$size))
      if(length(tryagain)>0){
        try(utils::download.file(url = addrs[tryagain], destfile = tmp.files[tryagain], method = "libcurl", quiet = T, mode = "wb"))
      }
      
      finfo<-file.info(tmp.files)
      tryagain<-which(finfo$size<30000000 | is.na(finfo$size))
      if(length(tryagain)>0){
        try(utils::download.file(url = addrs[tryagain], destfile = tmp.files[tryagain], method = "libcurl", quiet = T, mode = "wb"))
      }
      
    }
    
    
    if (isTRUE(ncores == 1)) { # Run in series
      if (verbose) cat("\nRunning in series...")
      `%how%` <- foreach::`%do%`
    }
    
    
  if (isTRUE(ncores > 1)) { # Run in parallel
    if (verbose) cat("\nRunning in parallel...")
    suppressWarnings(file.remove("BARRA_download_log.txt", showWarnings = FALSE))
    cl <- parallel::makeCluster(ncores, outfile = "BARRA_download_log.txt")
    doParallel::registerDoParallel(cl)
    if (verbose) cat(paste("\nProgress log output to:\n", getwd(), "/BARRA_download_log.txt", sep = ""))
    on.exit(expr = closeAllConnections())
    `%how%` <- foreach::`%dopar%`
  }

  all.yrs.weather <- foreach::foreach(y = seq_along(years), .combine = rbind, .multicombine = T, .export = "nc.process") %how% {
    if (verbose) cat("\nStarting", years[y],": ")
      all.mons.weather <- list()
      for (m in 1:length(mons)) {
        if (verbose) cat(month.abb[m], "|", sep = "")
        nc.path <- paste(tmp.dir,"/",paste(vars[v], years[y], mons[m], sep="_"), ".nc", sep = "")
        nc.data <- try(nc.process(nc.path))
        if(class(nc.data)=="try-error"){
          adrs<-paste("https://thredds.nci.org.au/thredds/fileServer/ob53/output/reanalysis/AUS-11/BOM/ERA5/historical/hres/BARRA-R2/v1/day/",
                      vars[v], "/latest/", vars[v], "_AUS-11_ERA5_historical_hres_BOM_BARRA-R2_v1_day_", years[y], mons[m], "-", years[y], mons[m], ".nc",
                      sep = ""
          )
          try(utils::download.file(url = adrs, destfile = nc.path, method = "libcurl", quiet = T, mode = "wb"))
          nc.data <- try(nc.process(nc.path))
        }
        file.remove(nc.path)
        all.mons.weather[[m]] <- nc.data
      }
      gc(full = T)
 
      tnames <- unlist(lapply(all.mons.weather, function(x) dimnames(x)[[3]]))
      all.mons.weather <- abind::abind(all.mons.weather, along = 3)

      env.info.yr.sub <- data.frame(
        "Environment" = Envs[Years == years[y]],
        "Lat" = Lats[Years == years[y]],
        "Lon" = Lons[Years == years[y]]
      )
      env.info.yr.sub$lon.ind <- sapply(env.info.yr.sub$Lon, function(x) which.min(abs(as.numeric(dimnames(all.mons.weather)[[1]]) - as.numeric(x))))
      env.info.yr.sub$lat.ind <- sapply(env.info.yr.sub$Lat, function(x) which.min(abs(as.numeric(dimnames(all.mons.weather)[[2]]) - as.numeric(x))))
      env.weather <- sapply(seq_len(nrow(env.info.yr.sub)), function(x) all.mons.weather[env.info.yr.sub$lon.ind[x], env.info.yr.sub$lat.ind[x], ])
      env.weather <- t(env.weather)
      #env.weather <- as.matrix(env.weather, row.names = env.info.yr.sub$Environment)[, 1:365]
      if (ncol(env.weather) < 365) {
        env.weather <- cbind(env.weather, matrix(NA, nrow = nrow(env.weather), ncol = 365 - ncol(env.weather)))
      }
      rownames(env.weather) <- env.info.yr.sub$Environment
      colnames(env.weather) <- 1:ncol(env.weather)
      return(env.weather)
  }
    if (isTRUE(ncores > 1)) { # Run in parallel
    parallel::stopCluster(cl)
    doParallel::stopImplicitCluster()

    if (verbose) {
      cat("\nFinished parallel :)\n")
    }
    Sys.sleep(2)
    suppressWarnings(file.remove("BARRA_download_log.txt"))
  }
  gc(full = T)
  
  all.yrs.weather<-all.yrs.weather[Envs,]
  if (verbose & sum(is.na(all.yrs.weather)) > 0) {
        NAenvs <- Envs[!complete.cases(all.yrs.weather)]
        cat("\nNAs returned at ", paste(NAenvs, collapse = " "))
  }
  
  all.vars.weather[[v]] <- all.yrs.weather
  }
  
  names(all.vars.weather) <- vars

  all.vars.weather$pr <- all.vars.weather$pr * 86400 # convert "kg m-2 s-1" to "mm day-1"
  all.vars.weather$tasmax <- all.vars.weather$tasmax - 273.15 # convert from Kelvin to deg C
  all.vars.weather$tasmin <- all.vars.weather$tasmin - 273.15 # convert from Kelvin to deg C
  all.vars.weather$rsdt <- all.vars.weather$rsdt / 41.67 * 3.6 # convert from w/m2 to MJ/m2

  names(all.vars.weather) <- c("daily_rain", "max_temp", "min_temp", "relhumidity", "radiation")


  # Calculate VPD
  all.vars.weather$vp_deficit <- vpdfun(
    tmin = all.vars.weather$min_temp,
    tmax = all.vars.weather$max_temp,
    relh = all.vars.weather$relhumidity
  )

  all.vars.weather <- all.vars.weather[!names(all.vars.weather) == "relhumidity"]

  DLs <- t(sapply(Lats, function(x) springpheno::daylength(daystop = 365,lat = x)))
  rownames(DLs) <- Envs
  all.vars.weather$day_length <- DLs

  env.info <- data.frame("Environment" = Envs, "Lat" = Lats, "Lon" = Lons)

  out <- list("data" = all.vars.weather, "Env.info" = env.info)
  return(out)
}
