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

###### Functions

download_observations(
  # loads up the Portal Data in your working folder.
  # this doesn't need to be run everytime unless you
  # are updating with newest data
  
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
  # used to find skipped censuses or two censuses
  # in one month.
  # requires data$month column in a date format
  # returns dataframe with time difference column
  
  timestep = diff(data$month)
  timestep = append(timestep, 40)
  data$mo_diff = timestep
  return(data)
}

cleanup_dataframe = function(data, type="rodent"){
  #  basic data frame clean up which differs between the rodent data and the
  # NDVI and weather data
  
  if(type=="rodent"){
    data = data |> filter(drop !="drop") |> select(!drop)
    data = data |> mutate(month = yearmonth(month))
  }
  data = data |> mutate(nmonth=month(month),year = year(month), time=row_number() )
  return(data)  
}

### rodent data functions
get_rodent_data = function(currency, lt_trt){
  # generates monthly cleaned data for all species on long-term plots (controls:4, 11, 14, 17  
  # or krat treatments:). Can specify abundance or energy use. 

  # summarize_rodent_data needs min_plots because plot 24 is missing in early years. 
  # download_if_missing will insert the missing newmoon months into the data frame.
  # even though there is no data so we can interpolate the missing info in later steps
  
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
  
  # for missing censuses, need to change NA for treatment to treatment data being 
  # generated.
  # first three censuses have more control plots so these are dropped
  # selects for the treatment passed to get_rodent_data()
  # adjusts counts by number of plots actually sampled
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
  # because of not-quite-monthly sampling. We need output to find gaps and 
  # multiple samples/month. This is so we can adjust the 'month' of sampling to
  # the month that the sampling event was intended to sample. Writes to files
  # where this is done by hand. This function is only needed if adj_yearmonths.csv 
  # does not already exist (i.e. this generates the precursor file that is then
  # adjusted to make adj_yearmonths.csv). Once adj_yearmonths.csv exists, there 
  # is a different process for updating that file by hand.

  # function requires data with columns: month, newmoonnumber, nplots censusdate
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

######## Data generating scripts

# Make control data
control_data = get_rodent_data("abundance", "control")
# write.csv(rodent_data, "raw_rodents.csv") # used for generate_adj_yearmonths_file.R

# Before running next step, update adj_yearmonth.csv.
# 1. Open control_data and compare newest dates with adj_yearmonth.csv
# 2. If newer dates in control_data, then add to adj_yearmonth.csv:
#  * newmoonnumber
#  * Census date: first night census date, 
#  * month: Year Mon(th) that census should be (i.e. was this the second trip in a month?
#       should this be the next month's census?)
#  * mo_diff: how many months since the previous census? 1 would be Aug then Sep,
#             2 would be Aug then Oct, 0 would be Aug then Aug.
#  * drop: sometimes two censuses do just occur in the same month. i.e. census pattern:
#             Aug Sep Sep Oct Nov... The second census isn't the Oct census because 
#             we trapped in Oct. Nor is the first one the August census, because we 
#             trapped then too. In that case, type "drop" to drop the second census 
#             in the month.

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
