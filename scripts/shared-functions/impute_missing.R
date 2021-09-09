# This function takes a data frame as its input, and uses mean imputation or linear interpolation depending on the data
impute_missing <- function(X = dat,
                           method = "min-impute",
                           replace.inf = T){
  
  # Find and replace Inf values with NA
  X <- X %>% 
    mutate(across(
      !c(iso3c),
      ~if_else(is.infinite(.x), as.numeric(NA), .x)
    ))
  
  #identify the variables with missing data
  missing <- X %>%
    select(
      where(
        ~any(is.na(.x))
      )
    ) %>%
    names()
  
  #identify which columns do not vary over time
  timeInvariant <- X %>%
    group_by(iso3c) %>%
    summarise(
      across(
        everything(),
        ~length(unique(.x))
      )
    ) %>%
    select(
      where(
        ~all(.x == 1)
      )
    ) %>%
    names()
  #the cols the change over time must be everything else minus iso3c and date
  timeVariant <- 
    setdiff(names(X), c(timeInvariant, "iso3c", "date"))
  
  #for the columns that do vary over time we use mean imputation weighted by
  #the distance calculated with the other time invariant variables
  mean_impute <- function(data, missing_vars, weight_vars, countryVar = "iso3c", timeVar = "date"){
    #calculate weights
    scaled <- data %>%
      select(all_of(countryVar), all_of(weight_vars)) %>%
      unique() %>%
      mutate(
        across(
          !all_of(countryVar),
          ~.x/sum(.x, na.rm = TRUE)
        )
      )
    #check which countries need imputing
    countries <- data %>% ungroup() %>%
      group_by(across(all_of(countryVar))) %>%
      summarise(
        across(
          all_of(missing_vars),
          ~any(is.na(.x))
        )
      ) %>%
      rowwise() %>%
      mutate(
        missing = any(c_across(missing_vars))
      ) %>%
      filter(
        missing
      ) %>% pull(all_of(countryVar))
    #set up data
    data <- data %>% mutate(distance = NA)
    if(timeVar %in% names(data)){
      data <- data %>% group_by(
        across(all_of(timeVar))
      )
    }
    #impute for a country
    for(country in countries){
      scaled_specific <- sweep( #subtract the country#s variables from every other variable
        as.matrix(
          select(scaled, !all_of(countryVar))
        ),
        2,
        as.numeric(
          select(filter(scaled, iso3c==country), !all_of(countryVar))
        )
      )^2 %>% #square, then sum across the rows
      rowSums(na.rm = TRUE)
      #now square root and invert
      scaled_specific <- (scaled_specific)^(-1/2)
      #set the value for our current country (will be inf) to 0
      scaled_specific[is.infinite(scaled_specific)] <- 0
      #add iso3c backin
      scaled_specific <- cbind(scaled_specific,
                               scaled[[countryVar]]) %>%
        as.data.frame()
      colnames(scaled_specific) <- c("distance", "iso3c")
      #impute values using these as weights
      data <- data %>%
        select(!distance) %>%
        left_join(scaled_specific, join = "iso3c") %>%
        mutate(
          distance = as.numeric(distance),
          across(
            all_of(missing_vars),
            ~if_else(
              is.na(.x) & iso3c == country,
              weighted.mean(.x, distance, na.rm = TRUE),
              .x
            )
          )
        )
    }
    return(data %>%
             ungroup() %>%
             select(!distance))
  }
  timeInvariant_missing <- intersect(missing, timeInvariant)
  X <- X %>%
    select(!all_of(timeInvariant_missing)) %>%
    left_join(
      X %>%
        select(all_of(timeInvariant), iso3c) %>%
        unique() %>%
        mean_impute(missing_vars = timeInvariant_missing, weight_vars = timeInvariant) %>%
        select(all_of(timeInvariant_missing), iso3c),
      by = "iso3c"
    )
  
  #for the time varying variables we use a more complex procedure
  #for daily rates/mobility we linearly interpolate to fill missing values, assuming rates
  #rise from 0 if the NA occurs before any data, and that rates remain the same 
  #for any NA that occur after the data we assume it stays at the last non-na value,
  #for cumulative or average values we calculate these once the other values are imputed
  #for countries with no values for a given variable we mean impute these at the end
  better_approx <- function(x, y, xout, rule){
    #a wrapper for approx to catch when there are no values at all
    if(sum(!is.na(y)) < 2){
      #if no values (or only one) just return them to be imputed else where
      return(
        y
      )
    } else{
      #before interpolating if a leading value is NA we set it to 0
      if(is.na(y[1])){
        y[1] <- 0
        #this means that our leading values are linear interpolated from 0 to whatever the first non NA value is
      }
      approx(x, y, xout, rule = rule)$y
      #by having rule = 2 any trailing NA get set to the last value
    }
  }
  timeVariant_missing <- intersect(missing, timeVariant)
  averages <- timeVariant_missing[str_detect(timeVariant_missing, "average")]
  cumulative <- timeVariant_missing[str_detect(timeVariant_missing, "cumulative")]
  lagged <- timeVariant_missing[str_detect(timeVariant_missing, "lagged")]
  vaccine_percentage <- timeVariant_missing[str_detect(timeVariant_missing, "vaccinated_pct")]
  timeVariant_missing <- setdiff(timeVariant_missing,
                                 c(averages, cumulative, lagged, vaccine_percentage))
  #we will linerly interpolate vaccination percentages, this doesn't make so much sense for keep the tailing NAs at the same value 
  #we cannot calculate these from vaccination rates due to the relation ship between 1/2 doses
  timeVariant_missing <- c(timeVariant_missing,
                             intersect(c("vaccinated_pct", "fully_vaccinated_pct"), vaccine_percentage))
  vaccine_percentage <- setdiff(vaccine_percentage, c("vaccinated_pct", "fully_vaccinated_pct"))
  #use linear interpolation on these values
  X <- X %>%
    group_by(iso3c) %>%
    arrange(iso3c, date) %>% #interpolate
    mutate(
      across(
        all_of(timeVariant_missing),
        ~better_approx(date,.x, date, rule = 2)
      )
    ) %>%
    ungroup()
  #check which still have NA
  timeVariant_missing <- X %>% select(all_of(timeVariant_missing)) %>%
    select(where(
      ~any(is.na(.x))
    )) %>% names()
  #mean impute on these values using the time invariant variables for weighting
  X <- X %>%
    select(!all_of(timeVariant_missing)) %>%
    left_join(
      X %>%
        select(all_of(timeVariant_missing), all_of(timeInvariant), iso3c, date) %>%
        mean_impute(missing_vars = timeVariant_missing, weight_vars = timeInvariant) %>%
        select(all_of(timeVariant_missing), iso3c, date),
      by = c("iso3c", "date")
    )
  
  #calculate the cumulative values
  #to keep these consistent and prioritise real data over imputed, we set any 
  #leading values to 0 but use the imputed values to fill the cumulative values forwards
  for(col in cumulative){
   X <- arrange(X, iso3c, date)
    #won't use across since its hard to access/use variable names
    daily_col <- str_remove(col, "cumulative_")
    #must do this by country
    countries <- X %>% 
      filter(is.na(.data[[col]])) %>% pull(iso3c) %>%unique()
    for(country in countries){
      if(all(is.na(X %>% filter(iso3c == country) %>% pull(all_of(col))))){
        #if no value just calculate from daily
        X[X$iso3c == country,col] <- cumsum(
          X[X$iso3c == country,daily_col]*(
            X[X$iso3c == country,"date"] - 
              lag(X[X$iso3c == country,"date"]$date, default = min(X[X$iso3c == country,"date"])-7)
            )
        )
      } else{
        if(is.na((X %>% filter(iso3c == country) %>% pull(all_of(col)))[1])){
          #fill leading zeros
          X[X$iso3c == country,][[col]][1:(which.min(is.na(X[X$iso3c == country,][[col]])) - 1)] <- 
            0
        }
        if(is.na(X %>% filter(iso3c == country) %>% pull(all_of(col)) %>% tail(1))){
          #use daily values to fill final points
          X[X$iso3c == country,][[col]][which.max(is.na(
            X[X$iso3c == country,][[col]]
          )):length(
            X[X$iso3c == country,][[col]]
          )] <- 
          X[X$iso3c == country,][[col]][which.max(is.na(
            X[X$iso3c == country,][[col]]
          )) - 1] + 
          cumsum(X[X$iso3c == country,][[daily_col]][which.max(is.na(
            X[X$iso3c == country,][[col]]
          )):length(
            X[X$iso3c == country,][[col]]
          )])
        }
      } 
    }
  }
  
  #calculate lagged values
  #most of these are just two leading values we we will set to 0, else we will
  #calculate from our data, won't quite match due to week etc but should be close enough
  first_weeks <- X %>%arrange(iso3c, date) %>% pull(date) %>% unique() %>% head(2)
  for(col in lagged){
    daily_col <- str_remove(col, "_lagged_two_weeks")
    X <- X %>%
      mutate(
        across(
          all_of(col),
          ~if_else(
            is.na(.x) & date %in% first_weeks,
            0,
            if_else(
              is.na(.x),
              lag(.data[[daily_col]], 2),
              .x
            )
          )
        )
      )
  }
  
  #calculate vacine percentage interactions
  #get rid of lagged variables
  vaccine_percentage <- vaccine_percentage[!str_detect(vaccine_percentage, "_lagged_two_weeks")]
  for(col in vaccine_percentage){
    #split into component variables
    vars <- str_split(col, "_over_")[[1]]
    if(vars[2] == "pop_65"){
      vars[2] <- "aged_65_older"
      X <- X %>% 
        mutate(across(
          all_of(col),
          ~if_else(is.na(.x),
                   (1-.data[[vars[1]]])/.data[[vars[2]]],
                   .x)
        ))
    } else{
      #three way interaction
      if(vars[2] == "covid_deaths"){
        vars[2] <- "daily_covid_deaths_per_100k"
      } else if(vars[2] == "covid_cases"){
        vars[2] <- "daily_covid_cases_per_100k"
      } else{
        stop(paste0(
          vars[2], " not recognised, please edit vaccine percentages section of impute_missing.R"
        ))
      }
      if(vars[1] == "vaccinated_pct_65_plus"){
        vars[1] <- "vaccinated_pct_over_pop_65"
      } else if(vars[1] == "fully_vaccinated_pct_65_plus"){
        vars[1] <- "fully_vaccinated_pct_over_pop_65"
      } else if(!vars[1] %in% c("vaccinated_pct", "fully_vaccinated_pct")){
        stop(paste0(
          vars[1], " not recognised, please edit vaccine percentages section of impute_missing.R"
        ))
      }
      X <- X %>% 
        mutate(across(
          all_of(col),
          ~if_else(is.na(.x),
                   (1-.data[[vars[1]]])*.data[[vars[2]]],
                   .x)
        ))
    }
  }
  
  #calculate averages, all but contiguous and excess death related we will calculate with weighted mean
  excess_averages <- averages[str_detect(averages, "excess")]
  averages <- setdiff(averages, excess_averages)
  sub_region_averages <- averages[str_detect(averages, "sub_region_average")]
  econ_region_averages <- averages[str_detect(averages, "econ_region_average")]
  averages <- setdiff(averages, c(sub_region_averages, econ_region_averages))
  region_averages <- averages[str_detect(averages, "region_average")]
  dist_averages <- averages[str_detect(averages, "dist_average")]
  contiguous_averages <- averages[str_detect(averages, "contiguous_country_average")]
  #calculate sub region averages
  for(var in sub_region_averages){
    main_var <- str_remove(var, "_sub_region_average")
    X <- X %>% 
      left_join(
        X %>% 
          group_by(date, 
                   across(
                     starts_with("subregion"))) %>%
          summarise(across(
            all_of(main_var),
            mean,
            .names = "mean_value"
          ))
      ) %>%
      mutate(
        across(
          all_of(var),
          ~if_else(
            is.na(.x),
            mean_value,
            .x
          )
        )
      ) %>%
      select(!mean_value)
  }
  #calculate econ region averages
  for(var in econ_region_averages){
    main_var <- str_remove(var, "_econ_region_average")
    X <- X %>% 
      left_join(
        X %>% 
          group_by(date, 
                   across(
                     starts_with("economist_region"))) %>%
          summarise(across(
            all_of(main_var),
            mean,
            .names = "mean_value"
          ))
      ) %>%
      mutate(
        across(
          all_of(var),
          ~if_else(
            is.na(.x),
            mean_value,
            .x
          )
        )
      ) %>%
      select(!mean_value)
  }
  #calculate region averages
  for(var in region_averages){
    main_var <- str_remove(var, "_region_average")
    X <- X %>% 
      left_join(
        X %>% 
          group_by(date, 
                   across(
                     starts_with("region"))) %>%
          summarise(across(
            all_of(main_var),
            mean,
            .names = "mean_value"
          ))
      ) %>%
      mutate(
        across(
          all_of(var),
          ~if_else(
            is.na(.x),
            mean_value,
            .x
          )
        )
      ) %>%
      select(!mean_value)
  }
  #calculate distance based averages
  for(var in dist_averages){
    main_var <- str_remove(var, "_dist_average")
    #which countries and times are missing
    missing_locations <- X %>%
      group_by(date, iso3c) %>%
      summarise(across(
        all_of(var),
        is.na
      )) %>% 
      filter(across(all_of(var))) %>% 
      arrange(iso3c, date) %>%
      select(date, iso3c)
    for(country in unique(missing_locations$iso3c)){
      #calculate weights based on longitude and lattiude
      country_loc <- X %>% 
        filter(iso3c == country) %>% 
        select(centroid_latitude, centroid_longitude) %>%
        unique()
      weights <- X %>%
        select(iso3c, centroid_latitude, centroid_longitude, is_subregionYES, population) %>%
        unique() %>%
        rowwise() %>%
        mutate(
          #use the haversine formula to calculate distance
          distance = pracma::haversine(
            c(centroid_latitude, centroid_longitude),
            c(country_loc$centroid_latitude, country_loc$centroid_longitude)
          ),
          weight = log(population)/log(distance+1),
          #if its a subregion we won't use it
          weight = if_else(
            is_subregionYES == 1 | iso3c == country,
            0,
            weight
          )
        ) %>%
        select(iso3c, weight)
      #calculate missing values
      X <- X %>% 
        group_by(date) %>%
        left_join(weights) %>%
        mutate(
          across(all_of(var),
          ~if_else(
            is.na(.x),
            weighted.mean(.data[[main_var]], weight),
            .x
            )
          )
        ) %>%
        select(!weight) %>%
        ungroup()
    }
  }
  #mean impute the remaining averages
  X <- X %>%
    select(!all_of(c(excess_averages, contiguous_averages))) %>%
    left_join(
      X %>%
        select(all_of(c(excess_averages, contiguous_averages)), all_of(timeInvariant), iso3c, date) %>%
        mean_impute(missing_vars = c(excess_averages, contiguous_averages), weight_vars = timeInvariant) %>%
        select(all_of(c(excess_averages, contiguous_averages)), iso3c, date),
      by = c("iso3c", "date")
    ) %>%#set any missing values for tests to 0
    mutate(
      across(
        all_of(c(excess_averages, contiguous_averages)),
        ~if_else(
          date == min(date) & is.na(.x),
          0,
          .x
        )
      )
    )
  #check what is left NA
  X %>%
    summarise(
      across(
        everything(),
        ~any(is.na(.x))
      )
    ) %>%
    select(where(~.x)) %>%
    names
  #for remaining excess death averages use linear interpolation to fill last values,
  #these will be places with no values at all
  #use linear interpolation on these values
  X <- X %>%
    group_by(iso3c) %>%
    arrange(iso3c, date) %>% #interpolate
    mutate(
      across(
        all_of(excess_averages),
        ~better_approx(date,.x, date, rule = 2)
      )
    ) %>%
    ungroup()
  #check what is left NA
  missing_still <-  
    X %>%
    summarise(
      across(
        everything(),
        ~any(is.na(.x))
      )
    ) %>%
    select(where(~.x)) %>%
    names()
  if(length(missing_still) != 0){
    warning(paste0(
      "The follow variables still have missing information: ",
      missing_still
    ))
  }
  return(X)
}
