#Prep data for calculations of first front, growing season, last frost, ave temperatures, etc. ----
prep1 <- function(data){
  if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse')
  if (!require('lubridate')) install.packages('lubridate'); library('lubridate')
  data$DATE <- ymd(data$DATE)
  data$dayofyear <- yday(data$DATE)
  data$YEAR <- year(data$DATE)
  data$MONTH <- month(data$DATE)
  data$DAY <- day(data$DATE)
  return(data)
}

#Calculate First Frost for one year at one site ----
firstfrost2 <- function(data) {
  data$truefalse <- ifelse(data$TMIN > 0, "TRUE", "FALSE") # For entire dataset, dd a column called truefalse with TRUE for positive temps and FALSE for negative temps
  data$truefalse[which(data$MONTH < 7)] <- NA # for all months before July, turn TRUE/FALSE into NA, so we only get fall freezing temps
  row <- suppressWarnings(min(which(data$truefalse == FALSE))) # return the row that corresponds to the minimum row of a freezing temp
  day <- data$dayofyear[row] #return the day of year for minimum row
  # day[is.infinite(day)] = NA # if day is infinite, replace with NA
  temp <- data$TMIN[row] # return the temperature for minimum row
  out <- c(day, temp)     # make a vector with the first frost day and the temperature on that day.
  return(out)
}

#Calculate First Frost for all years at one site ----
frostNEW <- function(data, site){
  listyear <- split(data, data$YEAR)
  out <- lapply(listyear, firstfrost2)
  df <- data.frame(matrix(unlist(out), nrow=length(out), byrow=T))
  colnames(df) <- c("firstfrostday", "firstfrosttemp")
  out2 <- cbind(site = rep(site, length(names(out))), 
                year = names(out), 
                df)
  return(out2)
}

#Daily minimum and maximum temperature throughout year ----
avetempsfn <- function(data, site){
  tempout <- data %>% group_by(dayofyear) %>%
    summarise(ave_min = mean(TMIN, na.rm = T),
              sd_min = sd(TMIN, na.rm = T),
              ave_max = mean(TMAX, na.rm = T),
              sd_max = sd(TMAX, na.rm = T))
  tempout2 <- cbind(site = rep(site, (nrow(tempout))), tempout)
  return(tempout2)
}

#Degree days ----


#Return daylength for a specific set of days in a year at one latitude ----
#NOTE: or you can just use the daylength function in the geosphere package!
daylengthfn <- function(latitude, days, site){
  j.constant = pi/182.625
  Axis.radians = 23.439 * pi/180
  
  daylength = 12 * (1 - tan(latitude * pi/180) * tan(Axis.radians * cos(j.constant * days)))
  daylength.corrected <- data.frame(day.old = seq(1,365), daylength = daylength) %>% 
    mutate(day.new = day.old - 10) %>%
    mutate(day.new = replace(day.new, day.new == -9, 356)) %>%
    mutate(day.new = replace(day.new, day.new == -8, 357)) %>%
    mutate(day.new = replace(day.new, day.new == -7, 358)) %>%
    mutate(day.new = replace(day.new, day.new == -6, 359)) %>%
    mutate(day.new = replace(day.new, day.new == -5, 360)) %>%
    mutate(day.new = replace(day.new, day.new == -4, 361)) %>%
    mutate(day.new = replace(day.new, day.new == -3, 362)) %>%
    mutate(day.new = replace(day.new, day.new == -2, 363)) %>%
    mutate(day.new = replace(day.new, day.new == -1, 364)) %>%
    mutate(day.new = replace(day.new, day.new == 0, 365))
  daylength.out <- daylength.corrected %>% select(c(daylength, day.new)) %>%
    relocate(day.new, .before = daylength) %>%
    arrange(day.new)
  daylength.out2 <- cbind(site = rep(site, (nrow(daylength.out))), daylength.out)
    return(daylength.out2)
}

