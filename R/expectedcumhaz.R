
# Expected cummulative hazard, age in days

expectedcumhaz <- function(ratetable, age, year, sex, time, method="exact", subdivisions = 100)
{
  if(method=="exact") {
    .year <- as.numeric(format( as.Date((1:time) + year, origin = "1960-01-01"), "%Y" ) ) 
    .temp <- diag(
      ratetable[as.character( pmin( floor( (age + (1:time) )/365.24), max(as.numeric(names(ratetable[, "2000", "male"]))) ) ),
            as.character( pmin( .year, max(as.numeric(names(ratetable["51", , "male"]))) ) ),
            sex] )
    
    return(sum(.temp))
  }
  
  if(method=="trapezoidal") {
    .f <- function(x)   { expectedhaz(ratetable=ratetable, age=age, year=year, sex=sex, time=x)}
    
    integrateA <- function(f, lower, upper, ..., subdivisions=100L, rel.tol=.Machine$double.eps^0.25,
                           abs.tol=rel.tol, stop.on.error=TRUE, keep.xy=FALSE, aux=NULL)
    {
      r <- integrate(f, lower, upper, ..., subdivisions=subdivisions, rel.tol=rel.tol,
                            abs.tol=abs.tol, stop.on.error=F, keep.xy=keep.xy, aux=aux)
      if ( !(r[['message']] %in% c('OK', 'maximum number of subdivisions reached')) ) {
        if (stop.on.error) { stop(r$message) }  }
      return(r)
    }
    
    return(integrateA(Vectorize(.f), lower=0, upper=time, subdivisions = subdivisions)$value)
  }
  
  if (method == "table") {
    
    table_scalar <- function(t) {
      birth_md <- format(
        as.Date(year - age, origin = "1960-01-01"),
        "%m-%d"
      )
      if (birth_md == "02-29") birth_md <- "03-01"
      
      bday <- as.Date(
        paste0(
          format(as.Date(year, origin = "1960-01-01"), "%Y"),
          "-", birth_md
        )
      )
      
      next_y <- as.Date(
        paste0(
          as.numeric(format(as.Date(year, origin = "1960-01-01"), "%Y")) + 1,
          "-01-01"
        )
      )
      
      end_date <- as.Date(year + t, origin = "1960-01-01")
      
      dates <- sort(unique(c(
        as.Date(year, origin = "1960-01-01"),
        bday + 365.24 * seq_len(200),
        next_y + 365.24 * seq_len(200),
        end_date
      )))
      dates <- dates[dates <= end_date]
      
      delta <- as.numeric(diff(dates))
      
      age_i  <- floor((age + cumsum(c(0, delta[-length(delta)]))) / 365.24)
      year_i <- as.numeric(format(dates[-length(dates)], "%Y"))
      
      age_i  <- pmin(age_i,  max(as.numeric(dimnames(ratetable)[[1]])))
      year_i <- pmin(year_i, max(as.numeric(dimnames(ratetable)[[2]])))
      
      haz <- ratetable[
        as.character(age_i),
        as.character(year_i),
        sex
      ]
      
      sum(haz * delta)
    }
    
    results <- vapply(time, table_scalar, numeric(1))
    
    return(matrix(results, ncol =1))
  }
  
}


