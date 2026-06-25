
# Expected instantaneous hazard, age in days
# date in number of days since 1960 / age in years / time in days

expectedhaz <- function(ratetable, age, year, sex, time, ...) 
{
  time <- pmin(time, 1000000)
  .year <- as.numeric(format( as.Date(time + year, origin = "1960-01-01"), "%Y" ) ) 
  .age <- round((age+time)/365.241) 
  idx <- cbind(as.character(.age), as.character(.year), sex
  )
  
  ratetable[idx]
}
  