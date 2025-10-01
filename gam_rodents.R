library("fpp3")
library("gratia")
library("mgcv")
library("tidyverse")
library("ggplot2")
library(patchwork)
#library("mgcViz")
library("DHARMa")
#library(itsadug)

model_checking = function(model_output){
  gam.check(model_output)
  acf(resid(model_output))
  pacf(resid(model_output))
  resids <- simulateResiduals(fittedModel = model_output, plot=TRUE)
}
# defining the seasonal knows used for all analyses

knots = list(c(.5,12.5))

# Rodent Data - counts of individuals (across all species)
data = read_csv("N_timeseries.csv")
data = data |> mutate(month=yearmonth(month))

# see the data
ggplot(data, aes(x=month,y=total)) +
  geom_line()
# Season and Year independent
total_s_y.m = gam(total ~ s(nmonth, k=12,bs="cc") + s(year, k=45, bs='cr'), 
              data=data, method="REML", family="tw", knots=knots)

# check residuals and smooths
model_checking(total_s_y.m)
summary(total_s_y.m)

# season and year interaction included
total_sxy.m = gam(total ~ s(nmonth, k=12,bs="cc") + s(year, k=40, bs='cr') + ti(nmonth, year, bs=c("cc","cr")), 
                 data=data, family="tw", method="REML", knots=knots)

model_checking(total_sxy.m)
summary(total_sxy.m)


# compare models
anova(total_s_y.m, total_sxy.m, test="F")

# viz best model
draw(total_sxy.m, scales="free", rug=FALSE, n=50) + plot_layout(widths =1) &
  theme(strip.text.x = element_text(size=8))

# NDVI - monthly greeness index, bounded by 0 and 1
NDVI = read_csv("NDVI.csv", col_types="ncdc")
NDVI = NDVI |> select(-c(1)) 
NDVI = NDVI |> mutate(month=yearmonth(month), nmonth=month(month), year = year(month))

# see data
ggplot(NDVI, aes(x=month,y=ndvi)) +
  geom_line()

# NDVI season and year independently
ndvi_s_y.m = gam(ndvi ~ s(nmonth, k=12,bs="cc") + s(year, k=40, bs='cr'), 
              data=NDVI, method="REML", family=inverse.gaussian(link = "log"), knots=knots)
model_checking(ndvi_s_y.m)
summary(ndvi_s_y.m)

# NDVI - season and year interaction included

ndvi_sxy.m = gam(ndvi ~ s(nmonth, k=12,bs="cc") + s(year, k=40, bs='cr') + ti(nmonth, year, bs=c("cc","cr")), 
                 data=NDVI, family=inverse.gaussian(link = "log"), method="REML", knots=knots)
model_checking(ndvi_sxy.m)
summary(ndvi_sxy.m)

# compare models
anova(ndvi_s_y.m, ndvi_sxy.m, test="F")

# viz best model
draw(ndvi_sxy.m, scales="free", rug=FALSE, n=50) + plot_layout(widths =1) &
  theme(strip.text.x = element_text(size=8))

# Weather rain and temperature - monthly rain in mm, temperature in C

PPT = read_csv("ppt_T.csv", col_types="nnnddddc")
PPT = PPT |> select(-c(1))

# Rain
ggplot(PPT, aes(x=yearmonth(month),y=precipitation)) +
  geom_line()

rain_s_y.m = gam(precipitation ~ s(month_col, k=12,bs="cc") + s(year, k=45, bs='cr'), 
             data=PPT, method="REML", family=tw(), knots=knots)
model_checking(rain_s_y.m)
summary(rain_s_y.m)

# season and year interaction included
rain_sxy.m = gam(precipitation ~ s(month_col, k=12,bs="cc") + s(year, k=45, bs='cr') + ti(month_col, year, bs=c("cc","cr")), 
                data=PPT, family="tw", method="REML", knots=knots)
model_checking(rain_sxy.m)
summary(rain_sxy.m)

anova(rain_s_y.m,rain_sxy.m, test="F")

# viz best model
draw(rain_s_y.m, scales="free", rug=FALSE, n=50) + plot_layout(widths =1) &
  theme(strip.text.x = element_text(size=8))

# Temperature

ggplot(PPT, aes(x=yearmonth(month),y=meantemp)) +
  geom_line()

# season and year independent
temp_s_y.m = gam(meantemp ~ s(month_col, k=12,bs="cc") + s(year, k=45, bs='cr'), 
            data=PPT, method="REML", family="gaussian", knots=knots)
model_checking(temp_s_y.m)
summary(temp_s_y.m)

# season and year interaction included
temp_sxy.m = gam(meantemp ~ s(month_col, k=12,bs="cc") + s(year, k=40, bs='cr') + ti(month_col, year, bs=c("cc","cr")), 
               data=PPT, family="gaussian", method="REML", knots=knots)
model_checking(temp_sxy.m)
summary(temp_sxy.m)

anova(temp_s_y.m, temp_sxy.m, test="F")

# viz best model
draw(temp_sxy.m, scales="free", rug=FALSE, n=50) + plot_layout(widths =1) &
  theme(strip.text.x = element_text(size=8))



# Story from the results - graphs are from the best fit model

# Precipitation 
# trend: precipitation declines from 1980-2010ish then plateaus
# seasonality: trough in may, peak in august, stable from oct-feb.
# interaction: none. No change in the seasonal signal of ppt over time. Maybe.
# So the interaction term is significant but model selection suggests not enough to warrant 
# the increase in df?

ggplot(PPT, aes(x=yearmonth(month),y=precipitation)) +
  geom_line()
draw(ppt.m, scales="free", rug=FALSE, n=50) + plot_layout(widths =1) &
  theme(strip.text.x = element_text(size=8))


# Temperature -
# trend: temperature shifts in mid 1980s, rises util 2000, then plateaus
# seasonality: peaks in june/july
# interaction: I believe the interaction is related to the change in 2000 but remains stable afterwards?

ggplot(PPT, aes(x=yearmonth(month),y=meantemp)) +
  geom_line()
draw(temp.full, scales="free", rug=FALSE, n=50) + plot_layout(widths =1) &
  theme(strip.text.x = element_text(size=8))

# NDVI (remember model is inverse): 
# Trend: NDVI stable until 2010, declines in 2010, new highs after 2015?
# Seasonality: highest in august
# Interaction: yes. If I understand this right, then 1990-2000 winter and summer greeneess 
# peaks were closer together. 2000-2020 summers increased productivity but winters diminished.
# potential signs that after 2020 we have begun to flip.

ggplot(NDVI, aes(x=month,y=ndvi)) +
  geom_line()
draw(ndvi.full, scales="free", rug=FALSE, n=50) + plot_layout(widths =1) &
  theme(strip.text.x = element_text(size=8))

# rodent abundance (remember model is inverse):
# trend: increasing, may have flattenedin 2000. Probably need a AR model to decrease wiggiligness
# seasonality: highest in june/july
# interaction: significant. 1980-2000 - less seasonal? After 2000 more seasonal.Also, potential shift in timing 
# from 1977-1990 - DS effect?

ggplot(data, aes(x=month,y=total)) +
  geom_line()
draw(total.full, scales="free", rug=FALSE, n=50) + plot_layout(widths =1) &
  theme(strip.text.x = element_text(size=8))

# Questions - why are the partials inverse from what I'm expecting????? 
# Answer: Something about gamma? I redid models with gaussian and the viz is literally the inverse of the gamma model


# Next steps:
  # 1) need to put in autocorrelation. Will help clarify the trends by removing autocorrelation
  # 2) some of the data models are pretty piss. Need to figure out better distribution families to use
  # 3) add in plant census data because plant story may be compositional - shrubs data and quadrat data.




