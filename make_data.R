#----------------------------
# Gets monthly Portal rodent, NDVI, rain, and temperature
#   data. Fills missing months. Rodent data from 
#   long-term controls
#---------------------------




library(lubridate)
library(portalr)
library(portalcasting)
library(fabletools)
library(tsibble)
library(tidyverse)

download_observations(
  path = ".",
  version = "latest",
  source = "github",
  quiet = FALSE,
  verbose = FALSE,
  pause = 30,
  timeout = getOption("timeout"),
  force = TRUE
)

calc_diff = function(data){
  # calculate time difference between rows
  # requires data$month column in a date format
  # returns dataframe with time difference column
  
  timestep = diff(data$month)
  timestep = append(timestep, 40)
  data$mo_diff = timestep
  return(data)
}

### RODENT DATA
get_rodent_data = function(currency, lt_trt){
  # extracts long-term control plots (4, 11, 14, 17) for every species, need min_plots because plot
  # 24 is missing in early years. This combo will give us the missing newmoon months.
  
  abundance_data <- portalr::summarize_rodent_data(path = ".", 
                                                   clean = FALSE,
                                                   level = "Treatment", 
                                                   plots = "Longterm",
                                                   type= "Rodents",
                                                   unknown = FALSE,
                                                   shape = "crosstab",
                                                   output = currency,
                                                   time = "all",
                                                   na_drop = FALSE,
                                                   zero_drop = FALSE,
                                                   min_traps=1,
                                                   effort = TRUE,
                                                   min_plots = 1,
                                                   download_if_missing = TRUE,
                                                   quiet=FALSE)
  
  # missing censuses need treatment added, first three censuses have more plots, adjusts counts
  # by number of plots sampled, formats date
    rodent_abundance = abundance_data |> 
    replace_na(list(treatment=lt_trt)) |>
    filter(newmoonnumber > 3)|> 
    filter(treatment == lt_trt) |> 
    mutate(month=yearmonth(censusdate)) |>
    mutate(across(BA:SO, ~.x*(4/nplots))) 
  
  # interpolate missing census values, calculate total abundance
  
   newmoon_abundance = rodent_abundance|> 
    mutate(across(BA:SO, portalcasting::round_na.interp)) |>
    rowwise() |> mutate(total = sum(c_across(BA:SO))) |>
    ungroup()
  return(newmoon_abundance)
}

make_yearmonth_gaps = function(data, filename){
  # needs output from get_rodent_data() to find gaps and multiple samples/month
  # in time series
  # requires data with columns: month, newmoonnumber, nplots censusdate
  # month column must be in yearmonth format
  
  match_yearmonths = calc_diff(data) |> select(newmoonnumber, nplots, censusdate, month, mo_diff)
  write.csv(match_yearmonths, filename)
}

merge_newmonths = function(data, filename){
  # Merges adjusted yearmonths with datafile and cleans up columns
  adj_yearmonth = read.csv(filename)
  month_abundance = left_join(data, adj_yearmonth, by='newmoonnumber')
  month_abundance = month_abundance |> select(!c(censusdate.y,month.x, mo_diff))
  month_abundance = month_abundance |> rename( censusdate = censusdate.x, month = month.y)
  return(month_abundance)
}


cleanup_dataframe = function(data, type="rodent"){
   if(type=="rodent"){
     data = data |> filter(drop !="drop") |> select(!drop)
   }
  data = data |> mutate(nmonth=month(month),year = year(month), time=row_number() )
return(data)  
}

# Make control data
control_data = get_rodent_data("abundance", "control")
# write.csv(rodent_data, "raw_rodents.csv") # used for generate_adj_yearmonths_file.R
monthly_data = merge_newmonths(control_data, "adj_yearmonths.csv")
output = cleanup_dataframe(monthly_data, "rodent")
write.csv(output, "N_timeseries_control.csv", row.names = FALSE) # output file for analysis

# Make control energy data
controlE_data = get_rodent_data("energy", "control")
# write.csv(rodent_data, "raw_rodents.csv") # used for generate_adj_yearmonths_file.R
monthlyE_data = merge_newmonths(controlE_data, "adj_yearmonths.csv")
output = cleanup_dataframe(monthlyE_data, "rodent")
write.csv(output, "E_timeseries_control.csv", row.names = FALSE) # output file for analysis


# Make krat exclosure data
exclosure_data = get_rodent_data("abundance", "exclosure")
# write.csv(rodent_data, "raw_rodents.csv") # used for generate_adj_yearmonths_file.R
monthly_data = merge_newmonths(exclosure_data, "adj_yearmonths.csv")
output = cleanup_dataframe(monthly_data, "rodent")
write.csv(output, "N_timeseries_exclosure.csv", row.names = FALSE) # output file for analysis
### NDVI DATA

get_ndvi_data = function(filename){
  # fetches NDVI data from Portal Database and processes for analysis
  # writes file. Filename is required input.
  
  month_ndvi = ndvi(
    level = "monthly",
    sensor = "landsat",
    fill = TRUE,
    path = ".",
    download_if_missing = TRUE
  )
  month_ndvi = month_ndvi |> mutate(month = yearmonth(date)) |>
    rename()
  write.csv(month_ndvi, filename)
}

get_ndvi_data("NDVI.csv")

### Precipitation data

get_weather_data = function(filename){
  month_weather = weather(
  level = "monthly",
  fill = TRUE,
  path = "."
  )
  month_weather = month_weather |> rename(month_col = month) |>
    select(year, month_col, meantemp, precipitation, anomaly_ppt, anomaly_meant) |> 
    mutate(month = yearmonth(paste(year,"-",month_col))) 
write.csv(month_weather,filename)
}


get_weather_data("ppt_T.csv")
