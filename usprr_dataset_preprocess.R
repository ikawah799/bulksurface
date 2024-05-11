# This program arranges input data for the bulk-surface snowmelt model (Ikawa et al., 2024, WRR) 

rm(list = ls())
# library required for data processing
library(tidyverse)
library(lubridate) # required for handling time data (e.g., make_datetime)
library(zoo)

######################################################################


apply_fun <- function(xnew, x, y) {
  predfunc <- approxfun(x[is.finite(x) & is.finite(y)],y[is.finite(x) & is.finite(y)], method="linear",rule=2)
  predfunc(xnew)
}

syear = 2011
eyear = 2020

T0 = 273.15
sb = 5.67*10^(-8) # Stefan-Boltzmann coefficient
ems = 0.95

setwd("C:/Users/Ikawa/Googledrive/Rproject/USPrr")

patheddy = "C:/Users/Ikawa/Googledrive/USPRR/w_qaqc_v2" # path for Ameriflux data
fileeddy = list.files(patheddy, full.names = T, pattern = ".csv") # fill.names = T adds directory of the file

dataall <- NULL

for (i in syear:eyear) {
  
  Year = as.integer(i)
  
 dataeddy0 <- as_tibble(read_csv(fileeddy[Year-2009]))
 dataeddy0[dataeddy0==-9999] = NA
 
 dataeddy1 <- dataeddy0 %>% 
 mutate(year  = floor(TIMESTAMP_END/1e8), 
       month = floor((TIMESTAMP_END-year*1e8)/1e6),
       day   = floor((TIMESTAMP_END-year*1e8 - month*1e6)/1e4),
       hour  = floor((TIMESTAMP_END-year*1e8 - month*1e6 - day*1e4)/1e2),
       mins  = floor((TIMESTAMP_END-year*1e8 - month*1e6 - 
                        day*1e4 - hour*1e2)),
       timestamp = make_datetime(year,month,day,hour,mins)) %>%
       mutate(doy = yday(timestamp - seconds(30)),
       year = year(timestamp-seconds(30)),
       month = month(timestamp-seconds(30)),
       hour = hour(timestamp-seconds(30))) %>%
       mutate(snowyear = if_else(doy >= 274, year+1, year)) %>% # DOY274 = Oct 1
       select(-day,  -mins)

 dataeddy2 <- dataeddy1 %>% 
 transmute(timestamp = timestamp, snowyear = snowyear, year = year, doy = doy, month = month, hour = hour, 
          ustar = USTAR_1_1_1, ws = WS_1_3_1, ws2 = WS_1_2_1, ws3 = WS_1_4_1, ta = TA_1_3_1, rh = RH_1_3_1, 
          rsd = SW_IN, rsu = SW_OUT, rld = LW_IN, rlu = LW_OUT, rld_floor = LW_BC_IN, rlu_floor = LW_BC_OUT, press = PA, prain = P_RAIN, 
          swc_1=SWC_1_1_1, swc_2 = SWC_1_2_1, swc_3 = SWC_1_3_1, swc_4 = SWC_1_4_1, swc_5 = SWC_1_5_1,
          ts_1 = TS_1_1_1, ts_2 = TS_1_2_1, ts_3 = TS_1_3_1, ts_4 = TS_1_4_1, ts_5 = TS_1_5_1, ts_6 = TS_1_6_1,
          ds_1 = D_SNOW_1_1_1, ds_2 = D_SNOW_1_1_2, ds_3 = D_SNOW_1_1_3,
          hflux = H_1_1_1+SH_1_1_1, leflux = LE_1_1_1+SLE_1_1_1, hflux_f = H_F_1_1_1, leflux_f = LE_F_1_1_1, hflux_f_floor = H_F_1_2_1, leflux_f_floor = LE_F_1_2_1, 
          gflux = (G_1_1_1 + G_1_1_2)/2, zL = ZL_1_1_1)
 
 dataall <- dataall %>% bind_rows(dataeddy2)
}

dataall2 <- dataall %>% 
  mutate(ws_f = if_else(is.finite(ws), ws, (ws2 + ws3)/2),
         tsr = ((rlu-(1-ems)*rld)/sb/ems)^(1/4)-T0,
         tsr_floor = ((rlu_floor-(1-ems)*rld_floor)/sb/ems)^(1/4)-T0,
         dsnow = (ds_1 + ds_2 + ds_3)/3,
         dsnow = ifelse(is.finite(dsnow),dsnow,(ds_2 + ds_3)/2))

# psnow calculation and fill small gaps
  dataall3 <- dataall2 %>% filter(snowyear >= syear & snowyear <= eyear) %>%
    filter((doy<=151) & (doy>=60)) %>% # March 1 to May 31
    group_by(snowyear) %>% group_split() %>%
  map_dfr(function(df){
    doy_complete <- 
      df %>% filter(dsnow < 1.) %>%
      slice_head(n = 1)
    df %>%
      mutate(
        SDD_dsnow = doy_complete$doy[1],
        dsnow = ifelse(doy > SDD_dsnow,0,dsnow),
        dsnow = na.approx(dsnow, rule=2),
        dsnow = ifelse(dsnow < 0,0,dsnow),
        dsnow_roll = rollapply(dsnow, 48, mean, align = "center", fill = NA),
        dsnow_roll = na.approx(dsnow_roll, rule=2, method = "constant"),
        psnow = c(0,diff(dsnow_roll)),
        psnow = if_else(psnow < 0, 0, psnow),
        
        ta_f = na.approx((ta), rule=2, maxgap=2),
        rh_f = na.approx((rh), rule=2, maxgap=2),
        rh_f = if_else(rh_f < 0, 0, rh_f),
        rh_f = if_else(rh_f > 100, 100, rh_f),
        ws_f = na.approx((ws_f), rule=2, maxgap=2),
        rsd_f = na.approx((rsd), rule=2, maxgap=2),
        rsu_f = na.approx((rsu), rule=2, maxgap=2),
        rld_f = na.approx((rld), rule=2, maxgap=2),
        rlu_f = na.approx((rlu), rule=2, maxgap=2),
        press_f = na.approx((press), rule=2, maxgap=2)
        )
  }) %>%
  select(-dsnow_roll) %>%
  mutate_if(is.numeric , replace_na, replace = NaN)

  # FURTHER GAPFILLING  
  dataall4 <- dataall3 %>% # Divide the period into 3 sections and fill by the mean 
    mutate(doycategory = if_else(doy<=90,1,if_else(doy<=120,2,3))) %>%
    group_by(snowyear, doycategory, hour) %>% group_split() %>%
    map_dfr(function(df){
      df %>%
        mutate(
          ta_f = if_else(is.nan(ta_f),mean(ta_f, na.rm=T),ta_f),
          rh_f = if_else(is.nan(rh_f),mean(rh_f, na.rm=T),rh_f),
          ws_f = if_else(is.nan(ws_f),mean(ws_f, na.rm=T),ws_f),
          rsd_f = if_else(is.nan(rsd_f),mean(rsd_f, na.rm=T),rsd_f),
          rsu_f = if_else(is.nan(rsu_f),mean(rsu_f, na.rm=T),rsu_f),
          rld_f = if_else(is.nan(rld_f),mean(rld_f, na.rm=T),rld_f),
          rlu_f = if_else(is.nan(rlu_f),mean(rlu_f, na.rm=T),rlu_f),
          press_f = if_else(is.nan(press_f),mean(press_f, na.rm=T),press_f),
        )
    }) %>%
    arrange(timestamp)
  
for (i in syear:eyear) {
  dataall_tmp <- dataall4 %>%
    select(timestamp, doy, SDD_dsnow, year, snowyear, month, hour, 
           hflux, leflux, hflux_f, leflux_f, hflux_f_floor, leflux_f_floor, gflux, ustar, zL, tsr, tsr_floor, dsnow, psnow,
           ta, rh, ws, rsd, rsu, rld, rlu, press, ta_f, rh_f, ws_f, rsd_f, rsu_f, rld_f, rlu_f, press_f,
           ts_1, ts_2, ts_3, ts_4, ts_5, ts_6) %>%
    mutate(flag = if_else(
      (snowyear == i), 1, 0)) %>% # March 1 to May 31
    filter(flag == 1) %>%
    select(-flag)
# readr::write_excel_csv(dataall_tmp, paste0("./snowin_usprr/snowinall",i,".csv"))
}
