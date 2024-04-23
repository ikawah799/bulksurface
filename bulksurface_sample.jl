# This program calculates energy budgets for a snow-covered forest (US-Prr & US-Uaf)
# according to the bulk surface approach
# Variables are explained in Ikawa et al (2024, Water Resources Research)

 using CSV, DataFrames, Plots, Dates, Statistics
 include("bulksurface_main.jl")
 
 targetyear = 2013
  
 rhosnowini = 170. # kg/m8
 ksoil_value = 0.8
 kmoss_value = 0.8
 dt = 30 * 60.

 zobs = 11.
 d = 1.46 # according to Chu et al (2021)
 z0 = 0.37 # 0.41  #2.0e-8
 zt = 3.34e-9 # 2.88e-6 #8.59e-7 #2.0e-8
 beta = 0.30 # 0.17
 dmoss = 0.14

 zobs_uaf = 8.
 z0_uaf = 0.65  #2.0e-8
 zt_uaf = 4.42e-5 # 5.02e-5 #2.0e-8
 beta_uaf = 0.22

 keco =     [1.5, 3.0, 1.9, 0.9, 0.7, 0.37, 0.6, 0.9,  0.47,  2.7]
 keco_uaf = [1.8, 2.5, 1.9, 1.2, 1.0, 0.5 , 1.65,  8.2, 1.52, 8.0]

# Copy to usprr_modelparams.R, usprr_modelerror.R
 SDD     = [119, 110, 133, 112, 110, 98, 118, 127, 100, 125]
 SDD_uaf = [113, 108, 131, 110, 105, 98, 115, 118,  93, 113]

  # observed soil temperature depth
 zsoilobs = [0.05, 0.1, 0.2, 0.3, 0.4, 1.0] .+ dmoss
 zsoilobs_uaf = [0.15, 0.3, 0.45, 0.8, 1.25]

 # bottom boundary of the soil simulation grid
 zbsoil = [0., .1, .3, .5, .7, .9, 1.1] .+ dmoss # do not change by site

 fidin = "bulksurface2013.csv"
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
tsoilini = soilinterp(zsoilobs, tgini, zbsoil, 7)
ksoil = zeros(length(zbsoil)) .+ ksoil_value
ksoil[1] = kmoss_value     
pswe = psnow/100 * rhosnowini/1e3
sweini = dsnow[1]/100 * rhosnowini/1e3

     (hflux_c, leflux_c, gflux_c, gflux0, Mf, teco_c, tgtotal_c, swetotal_c, dstotal_c, albedo_c, rhosnow, flag_snowmin, ustar_c)=
     bulksurface(ta, rh, ws, press, rsd, rld, psnow, zobs, d, z0, zt, beta, dt, ksoil,
     keco[Int64(targetyear-2010)], rhosnowini, sweini, tsoilini, zbsoil)

    plot(dstotal_c)