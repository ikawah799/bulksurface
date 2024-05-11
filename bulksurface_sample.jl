# This program calculates energy budgets for a snow-covered forest (US-Prr & US-Uaf)
# according to the bulk surface approach
# Variables are explained in Ikawa et al (2024, Water Resources Research)

# Output
# hflux: sensible heat flux (W m-2)
# leflux: latent heat flux (W m-2)
# gflux: heat flux from the bulk-surface to the snow layer (W m-2)
# gflux0: heat flux from the snow layer to the soil layer (W m-2)
# Mf: snowmelt rate (m s-1)
# teco: ecosystem surface temperature (deg C)
# tg: soil temperature (deg C)
# swe: snow water equivalent (m)
# dsnow: snow depth (m)
# albedo: albedo
# rhosnow: density of snow (kg m-3)
# flag_snowmin: flag to identify snow disapperance

# Input
# ta: air temperature (deg C)
# rh: relative humidity (%)
# ws: wind speed (m s-1)
# press: atmospheric pressrue (kPa)
# rsd: downward shortwave flux (W m-2)
# rld: downward longwave flux (W m-2)
# psnow: snow precipitation on the floor (m s-1 * time interval of input data)
# zobs: observation height for wind speed and air temperature (m) *if these are different, use different zm in "surfaceEB"
# d: displacement height (m)
# z0: roughness length for momentum (m)
# zt: roughness length for the heat (m)
# beta: evaporation efficiency
# dt: time step (s)
# ksoil: heat conductivity of soil (W m-1 K-1)
# keco: bulk ecosystem surface conductivity (W m-1 K-1)
# rhosnow_ini: initial snow density (kg m-3)
# swe_ini: initial snow water equivalent (m)
# tsoil_ini: initial soil temperature profile (deg C)
# zbsoil: bottom boundary depth for soil temperature (m)

 using CSV, DataFrames, Plots
 include("bulksurface_main.jl")
 
 targetyear = 2013
  
 rhosnow_ini = 170. # kg/m3
 ksoil_value = 0.8
 kmoss_value = 0.8
 dt = 30 * 60.

 zobs = 11.
 d = 1.46 # according to Chu et al (2021)
 z0 = 0.37
 zt = 3.34e-9
 beta = 0.30
 dmoss = 0.14

 zobs_uaf = 8.
 z0_uaf = 0.65
 zt_uaf = 4.42e-5
 beta_uaf = 0.22

# Bulk ecosysytem surface conductivity
 keco =     [1.5, 3.0, 1.9, 0.9, 0.7, 0.37, 0.6, 0.9,  0.47,  2.7]
 keco_uaf = [1.8, 2.5, 1.9, 1.2, 1.0, 0.5 , 1.65,  8.2, 1.52, 8.0]

# Snow disapperance date (SDD) based on albedo data
 SDD     = [119, 110, 133, 112, 110, 98, 118, 127, 100, 125]
 SDD_uaf = [113, 108, 131, 110, 105, 98, 115, 118,  93, 113]

  # observed soil temperature depth
 zsoilobs = [0.05, 0.1, 0.2, 0.3, 0.4, 1.0] .+ dmoss
 zsoilobs_uaf = [0.15, 0.3, 0.45, 0.8, 1.25]

 # bottom boundary of the soil simulation grid
 zbsoil = [0., .1, .3, .5, .7, .9, 1.1] .+ dmoss # same for both sites

 fidin = "usprr2013.csv"
 data=CSV.File(fidin,header=1,skipto=2,delim=",",missingstring="NA") |> DataFrame

ustar = data."ustar"
ws = data."ws_f"
ta = data."ta_f"
rh = data."rh_f"
rsd = data."rsd_f"
rsu = data."rsu_f"
rld = data."rld_f"
press = data."press_f"
psnow = data."psnow"
dsnow = data."dsnow"

ts = zeros(length(ta), 6)
ts[:, 1] = collect(Missings.replace(data."ts_1", NaN))
ts[:, 2] = collect(Missings.replace(data."ts_2", NaN))
ts[:, 3] = collect(Missings.replace(data."ts_3", NaN))
ts[:, 4] = collect(Missings.replace(data."ts_4", NaN))
ts[:, 5] = collect(Missings.replace(data."ts_5", NaN))
ts[:, 6] = collect(Missings.replace(data."ts_6", NaN))

tgini = ts[1,:]
tsoil_ini = soilinterp(zsoilobs, tgini, zbsoil, 7)
ksoil = zeros(length(zbsoil)) .+ ksoil_value
ksoil[1] = kmoss_value     
pswe = psnow/100 * rhosnow_ini/1e3
swe_ini = dsnow[1]/100 * rhosnow_ini/1e3

     (hflux, leflux, gflux, gflux0, Mf, teco, tg, swe, dsnow_model, albedo, rhosnow, flag_snowmin)=
     bulksurface(ta, rh, ws, press, rsd, rld, psnow, zobs, d, z0, zt, beta, dt, ksoil,
     keco[Int64(targetyear-2010)], rhosnow_ini, swe_ini, tsoil_ini, zbsoil)

# plot snow depth in cm
    plot(dsnow)
    plot!(dsnow_model*100)