# Bulksurface snowmelt model Ikawa et al (2024, WRR)
# Hiroki Ikawa (4/23/2024) 

using Interpolations
using Interpolations, GLM # used for initial soil temp profile

function bulksurface(ta::Vector{Float64}, rh::Vector{Float64}, ws::Vector{Float64}, 
    press::Vector{Float64}, rsd::Vector{Float64}, rld::Vector{Float64}, psnow::Vector{Float64}, 
     zobs::Float64, d::Float64, z0::Float64, zt::Float64, beta::Float64, dt::Float64, 
    ksoil::Vector{Float64}, keco::Float64, rhosnow::Float64, sweini::Float64, 
    tsoilini::Vector{Float64}, zbsoil::Vector{Float64}, site="prr"
)

# snow parameter
  lfusion = 333.6 * 1e3 # J/kg
  swemin = 1e-3 # 1mm based on swe

  rhosoil = 1.1e3
  rhowater = 1e3 # kg/m3
  rhoice = 920. # kg/m3
  rhomineral = 2650.
  rhoair = 1.29
 
  cwater = 4180.
  cice = 2100. # specific heat of snow J/K/kg
  cmineral = 870.
  cair = 1010.
  vwc = 1. - rhosoil/rhomineral
  crhosoil = cmineral * rhosoil + vwc* cice * rhoice +
  (1. - vwc - rhosoil/rhomineral) * cair * rhoair

  # bulk parameter
  ems = 0.95
  ldata = length(ta)
  nlayer = length(tsoilini) + 1
  flag_snowmin = Int64(1e10)
  flag_Mf = 0

  amax = 0.70
  amin = 0.30
  gamma_albedo = 5.19 # 0.837

  #if site=="uaf"
  #amax = 0.61
  #amin = 0.15
  #gamma_albedo = 4.6 # 0.837
  #end

    # OUTPUT 
    eflux = zeros(Float64, ldata)
    leflux = zeros(Float64, ldata)
    hflux = zeros(Float64, ldata)
    teco = zeros(Float64, ldata)
    tgtotal = zeros(Float64, ldata, nlayer)
    gflux = zeros(Float64, ldata)
    ustar = zeros(Float64, ldata)
    gflux0 = zeros(Float64, ldata)
    Mf = zeros(Float64, ldata) 
    swetotal = zeros(Float64, ldata)
    dstotal = zeros(Float64, ldata)
    albedo = zeros(Float64, ldata)
    rhosnowtotal = zeros(Float64, ldata)
    # Intermediate                                                                        
    zm = zeros(Float64, nlayer)
    ks = zeros(Float64, nlayer)
    crho = zeros(Float64, nlayer)

    nlayer = 8
    zbini = vcat(0, sweini[1]*rhowater/rhosnow, sweini[1]*rhowater/rhosnow .+ zbsoil)
    zm = zeros(Float64, nlayer)
    for k = 1:nlayer
        zm[k] = (zbini[k] + zbini[k+1]) * 0.5
    end

    tg =  vcat(ta[1],tsoilini)
    swe = copy(sweini[1])
    tgb = tg[end]

    count_snow = 0
    tgba = 0.

    for i = 1:ldata
        #
     # Input snowfall   
        swe = swe + psnow[i]/100 * rhosnow/rhowater

     # Determine albedo and input radiation rnplus
        
        albedo[i] = max(
            amin * amax * exp(gamma_albedo*swe*rhowater/rhosnow) /
                       (amax - amin + amin*exp(gamma_albedo*swe*rhowater/rhosnow))
                       ,0.3)

        if site == "uaf"
            albedo[i] = max(0.24 + 0.47 * swe*rhowater/rhosnow, amin)                        
        end

     rnplus = (1. - albedo[i]) * rsd[i] + ems*rld[i]
#    rnplus = rsd[i] - rsu[i] + ems*rld[i]
                 
     # Solve surface energy balance 
        (eflux[i], leflux[i], hflux[i], teco[i], gflux[i], ustar[i]) = surfaceEB(
            ta[i], rh[i], ws[i], press[i], rnplus, zobs,
            d, z0, zt, beta, tg[1], keco, zm[1], ems, flag_Mf
        )

       # if swe == swemin gflux[i] = 0. end
 
     # Soil scheme for i>1. i=1 is the boundary condition (ie, t=0)
        if i > 1           

             # Set the boundary and distance in between (zb, zm, dzm, dzb)
            zb = vcat(0., swe*rhowater/rhosnow, zbsoil .+ swe*rhowater/rhosnow)

            zm = zeros(nlayer)

            for k = 1:nlayer
                zm[k] = (zb[k] + zb[k+1]) * 0.5
            end

            dzm = diff(zb)
            dzb = cat(zm[1], diff(zm), dims = 1)
            ksnow = 10. ^ (0.00102 * rhosnow - 1.14) 
            ks[1] = ksnow
            ks[2:end] .= ksoil
            crho[1] = cice * rhosnow + (1. - rhosnow/rhoice) * cair * rhoair
            crho[2:end] .= crhosoil

            # ks[2] = 0.1 #*ksoil
            # crho[2] = 1.0*crhosoil

            # Specific heat times density for each layer
            for k = 2:nlayer
                if tg[k]<0. && tg[k]>-1.
                crho[k]=crho[k] + lfusion * vwc * rhowater
                end
            end
            
            # It would be more appropriate to provide teco rather than gflux here.
            # The fixed boundary condition is more stable than the flux boundary condition.
            # Solve for snow&soil temperature without considering snowmelt
            tg, = solvesoilkeco!(dzm, tg, teco[i], tgb, dt, crho, ks, keco, gflux[i], gflux[i-1])

            # Calculate snowmelt Mf
            if tg[1] > 0.
                dTpos = tg[1]
                tg[1] = 0.
            else
                dTpos = 0.
            end

            k0 = (1/ksnow*dzm[1]/2+1/ksoil[1]*dzm[2]/2)^(-1)*dzb[2]
            gflux0[i] = k0*(tg[1]-tg[2])/dzb[2]

            Mf[i] = dTpos * swe * cice / lfusion / dt # It is Mf devided by rhowater (m s-1)
            # Remind that "rhosnow * dsnow = rhowater * swe" in the snowmelt paper

            # Calculate the change in snow depth by fusion and sublimation
             if swe <= swemin
                dswe_fusion = 0.
                dswe_et = 0.
             else                
                dswe_fusion = Mf[i] * dt        # change in snow depth in m
                dswe_et = eflux[i] * dt / rhowater # change in snow depth in m
             end

             if dswe_fusion > 0.
                flag_Mf = 1
                swe = swe - dswe_fusion # max(dswe_fusion, dswe_et)
             else
                flag_Mf = 0
                swe = swe - dswe_et
             end

            if swe <= swemin
               swe = copy(swemin)
               flag_snowmin = min(flag_snowmin, i)
               break
            end
            
            # Bottom boundary scheme using Force-Restore model
            count_snow = count_snow + 1
            tgba = tgba + tg[end-1]

            if count_snow == Int64(24*3600/dt)
                tgba = tgba / (24*3600/dt)
                Dthermal = (1/ks[end-1] * dzm[end-1]/2 + 1/ks[end] * dzm[end]/2)^(-1) * dzb[end] / crho[end]
#                Dthermal = Dthermal * 1/4 # apply smaller thermal conductivity for the bottom (US-Prr specific)
                tgb = FRextended!(tgb, tgba, Dthermal, dzb[end],-1.)
                count_snow = 0
                tgba = 0.
            end
        end

        tgtotal[i, :] = tg
        swetotal[i, 1] = swe
        dstotal[i, 1] = swe.*rhowater./rhosnow

        rhosnowtotal[i] = rhosnow

        # rhosnow for the next time step
        rhosnow = rhosnow + rhosnow * swe/2.0 * rhowater * 9.8 /
        (0.3*1e-5*exp(0.02*rhosnow + 8110/(tg[1] + 273.15))) * dt
    
    end

    return hflux, leflux, gflux, gflux0, Mf, teco, tgtotal, swetotal, dstotal, albedo, rhosnowtotal, flag_snowmin, ustar
end

# Bulk surface energy balance scheme
# Hiroki Ikawa (10/21/2020) 

######################################################
# INPUTS
# ta     : air temp at zm (degee C)
# hm     : relative humidity at zm (#)
# ws     : wind speed at zm (m s-1)
# ustarm : measured ustar to constrain stability correction
# press  : pressure (k Pa)
# rnplus : net radiation + upward radiation (W m-2)
# zm     : observation height
# d      : displacement height
# z0     : aerodynamic roughness length
# zt     : roughness length for heat
# beta   : evaporation efficiency (=1 at water surface)
# t1     : temperature below the surface
# z1     : the distance between the surface and z1
# ks   : the ecosystem conductivity (W K-1 m-1)

# OUTPUTS
# eflux  : evapotranspiration (mm hr-1) 
# leflux : latent heat flux (W m-2)
# hflux  : sensible heat flux (W m-2)
# ts     : surface temperature (degree C)
# gflux  : below surface heat flux (W m-2)
# deltaT : ts - ta (degree C)
# count  : # of count required for interation 
# chu    : transfer velocity for sensible heat flux (m s-1)
# q      : specific humidity (g kg-1)
# ustar  : estimated ustar
######################################################

## test imput ###
#=
ta=25
t1=20
hm=70
ws=2
press=100 
rnplus=500
ustar=ws/10 
zm=11
d=2.91
zt=0.01
z0=0.2
beta=0.5
dz=0.01
ks=2.0

z0=1e-2
zt=1e-3
zq=1e-5
ems = 0.95
=#

    # LOWE (1976, Journal of Applied Meteorology) for eslowe (esat(T)) and deslowe (des/dT)
    eslowe(x) =
            6.107799961 +
            x * (
                4.436518521 * (10^-1) +
                x * (
                    1.428945805 * (10^-2) +
                    x * (
                        2.650648471 * (10^-4) +
                        x * (
                            3.031240396 * (10^-6) +
                            x * (2.034080948 * (10^-8) + x * 6.136820929 * (10^-11))
                        )
                    )
                )
            )
    deslowe(x) =
            4.438099984 * (1e-1) +
            x * (
                2.857002636 * (1e-2) +
                x * (
                    7.938054040 * (1e-4) +
                    x * (
                        1.215215065 * (1e-5) +
                        x * (
                            1.036561403 * (1e-7) +
                            x * (3.532421810 * (1e-10) - x * 7.090244804 * (1e-13))
                        )
                    )
                )
            )

function surfaceEB(ta::Float64, hm::Float64, ws::Float64, press::Float64, 
    rnplus::Float64, zm::Float64, d::Float64, z0::Float64, zt::Float64, beta::Float64, 
    t1::Float64, ks::Float64, dz::Float64, ems::Float64, flag_Mf::Int64)

    local  lv, ceu , chu , psim , eb    

    opt_steffensen = 1

    flag = 0
    # initial value
    ts = ta + 1e-2
    count = 0 # # of iteration
    cmax = 100 # maximum iteration
    hflux = 1e-10

    # constant and model setting
    sb = 5.67 * 10^(-8)  # Stefan-Boltzmann coefficient (W m-2 K-4)
    cp = 1.01 * 10^3     # specific heat of air (Jkg-1K-1)
    ck = 0.4           # von Karman constant
    grav = 9.8
    T0 = 273.15
    lfusion = 333.6 * 1e3 # J/kg

    srtcmu = ck / (log((zm - d) / z0))
    ustar = srtcmu * ws

    eps = 1e-4 # any small number. may be fine with 0.01 but not tested. 

    esat = eslowe(ta) # 6.1078 * exp(17.2694 * ta / (ta + 237.3)) # hPa
    e = esat * hm / 100.                       # vapor pressure hPa
    qa = 622. * e / (press * 10. - 0.378 * e) * 1e-3          # specific humidity kg/kg
    qasat = 622. * esat / (press * 10. - 0.378 * esat) * 1e-3          # specific humidity kg/kg
    dqT = deslowe(ta) * 0.622 * press * 10. / (press * 10. - 0.378 * e) .^ 2
    # dq/dT = (0.622*dqT*(ps - ...)/(ps - 0.378)^2

    rho  =
        1.293 * 273.15 / (273.15 + ta) *
        (press * 10. / 1013.25) *
        (1. - 0.378 * e / (press * 10.)) # air density  (about 1.2 kgm-3)  
    cprho = cp * rho
    Tv = (ta + T0) * (1. + 0.61 * qa) # virtual temperature 

    while true

        count = count + 1

        zoL = -ck * zm * grav * hflux / cprho / Tv / ustar^3

        if (zoL > 0.0001) && (zoL < 10.) && (count < cmax * 0.9)      # stable
            psih = 400. * log(1. + 7. / 400. * zoL + 0.005 * zoL^2)
            psim = 7. / 3. * log(1. + 3. * zoL + 10. * zoL^3)
        elseif (zoL < -0.0001) && (zoL > -10.) && (count < cmax * 0.9)      # unstable 
            xxh = sqrt(1. - 16. * zoL)
            psih = 2. * log(2. / (1. + xxh))
            xxm = sqrt(xxh)
            psim = log(8. / (xxm^2. + 1.) / (xxm + 1.)^2.) + 2. * atan(xxm) - 2. * atan(1)
        else
            psih = 0.
            psim = 0.
        end

        chu = ws * ck / (log((zm - d) / z0) + psim) * ck / (log((zm - d) / zt) + psih)
        if chu < 0.
            chu = ws * ck / (log((zm - d) / z0)) * ck / (log((zm - d) / zt))
        end

        ceu = chu*beta

       # ceu = ws * ck / (log((zm - d) / z0) + psim) * ck / (log((zm - d) / zq) + psih)
       # if ceu < 0.
       #     ceu = ws* ck / (log((zm - d) / z0)) * ck / (log((zm - d) / zq))
       # end

        if flag_Mf == 1 # snowmelt occurs
            lv = 2.5 * 1e6 - 2400. * ta       # latent heat for vaporization Jkg-1 (about 2.5 X 10^6 Jkg-1)  
        else flag_Mf == 0 # if snowmelt does not occur
            lv = 2.5 * 1e6 - 2400. * ta + lfusion # latent heat for vaporization Jkg-1 (about 2.8 X 10^6 Jkg-1)  
        end

        if opt_steffensen == 1
           # Steffensen method for the convergence 
           es = eslowe(ts) # 6.1078 * exp(17.2694 * ts / (ts + 237.3)) #hPa
           qs = 622. * es / (press * 10. - 0.378 * es) * 1e-3 # specific humidity kg/kg
           x1 =
              (rnplus - ks * (ts - t1) / dz - ems * sb * (ts + T0)^4. -
                  lv * rho * ceu * (qs - qa)) / (cprho * chu) + ta

          if x1 - ts > 1e2
              x1 = ts + 1e2
              flag = 1
          elseif x1 - ts < -1e2
              x1 = ts - 1e2
              flag = 1
          else
              flag = 0
          end

          es = eslowe(x1) # 6.1078 * exp(17.2694 * x1 / (x1 + 237.3)) #hPa
          qs = 622. * es / (press * 10. - 0.378 * es) * 1e-3 # specific humidity kg/kg

          x2 =
          (rnplus - ks * (x1 - t1) / dz - ems * sb * (x1 + T0)^4 -
              lv * rho * ceu * (qs - qa)) / (cprho * chu) + ta

          if abs(x2 - 2. * x1 + ts) < eps
              x3 = x2
              ts = x3
              break
          else
              x3 = ts - ((x1 - ts)^2.) / (x2 - 2. * x1 + ts)
              ts = x3
          end

          if count > cmax
              # fprintf('break')    
              break
          end
        else
            ts =
            (
                rnplus - lv * rho * ceu * (qasat - qa) - ems * sb * (ta + T0)^4. - ks * (ta - t1) / dz) /
            (cprho * chu + lv * rho * ceu * dqT + 4. * ems * sb * (ta + T0) ^ 3. + ks / dz
            ) + ta
            if count > 2
                break
            end
        end

        hflux = cprho * chu * (ts - ta) # required for zoL
    end

    es = eslowe(ts)  #hPa
    qs = 622. * es / (press * 10. - 0.378 * es) * 1e-3 # specific humidity g/kg
    leflux = lv * rho * ceu * (qs - qa) #Wm-2
    eflux = leflux / lv # * 3600 #mm hr-1 or kgH2O/m2/hr

    gflux = ks * (ts - t1) / dz
    deltaT = ts - ta

    srtcmu = ck / (log((zm - d) / z0) + psim)
    ustar = srtcmu * ws

    eb = (rnplus - ems * sb * (ts + T0)^4.) - (hflux + leflux + gflux)
 
    if abs(eb) > 100.
        println("too large eb")
     gflux = rnplus - ems * sb * (ts + T0)^4. - (hflux + leflux)
    end

    return eflux, leflux, hflux, ts, gflux, ustar, deltaT, count, chu, qa, eb, count

end

## Solve soil temperature profiles 
# Completely backward scheme: The balanced term in the prognostic equations is unknown.   
# The measurement location is assymetric within each grid, assuming the
# effect is small. 

# INPUT 
# zm     (m)        depth of the measurement
# tg    (deg C)     soil temperature 
# tg0   (deg C)     surface temperature 
# dt    (second)    time interval 
# ks    (W K-1 m-1) thermal conductivity 

# OUTPUT 
# tg    (deg C)     soil temperature 
# gflux (W/m2)      ground heat flux

# Calculated Variables
# zb    (m)         depth at grid boundary 
# dzb   (m)         delta z between zm(k) and zm(k-1)
# dzm   (m)         delta z at zm(k) 
# kb    (m2/s)      thermal diffusivity at grid boundary 

#---------------------------------------------------------
#         zb(k-1)  dzb(k-1)                   Ksb(k-1)
# zm(k-1)                   dzm(k-1)   Ks(k-1)
#         zb(k)    dzb(k)                     Ksb(k)
# zm(k)                     dzm(k)     Ks(k)
#         zb(k+1)  dzb(k+1)                   Ksb(k+1)
# zm(k+1)                   dzm(k+1)   Ks(k+1)
#---------------------------------------------------------

function solvesoilkeco!(dzm::Vector{Float64}, tg::Vector{Float64}, tg0::Float64, tgb::Float64, dt::Float64, 
    crho::Vector{Float64}, ks::Vector{Float64}, keco::Float64, gflux::Float64, gflux0::Float64)

    km = length(tg)
    zm = zeros(km)
    kb = zeros(km)

    zb = cat(0., cumsum(dzm), dims = 1)

    for k = 1:km
        zm[k] = (zb[k] + zb[k+1]) * 0.5
    end

    dzb = cat(zm[1], diff(zm), dims = 1)

    kb[1] = keco #Ks(1)
    for k = 2:km
#        kb[k] = (ks[k] * dzm[k]/2 + ks[k-1] * dzm[k-1]/2) /dzb[k]
     kb[k] = (1/ks[k] * dzm[k]/2 + 1/ks[k-1] * dzm[k-1]/2)^(-1) * dzb[k]
    end

    a = zeros(km)
    b = zeros(km)
    c = zeros(km)
    d = zeros(km)

    cn = 0.6

    for k = 2:km-1
        a[k] = cn * kb[k] / dzb[k]
        c[k] = cn * kb[k+1] / dzb[k+1]
        b[k] = a[k] + c[k] + dzm[k] * crho[k] / dt
        d[k] = tg[k] * dzm[k] * crho[k] / dt +
            (1 - cn) *
            (kb[k+1] * (tg[k+1] - tg[k]) / dzb[k+1] -
                kb[k] * (tg[k] - tg[k-1]) / dzb[k])
    end

    # The last two terms are flux divergent energy to melt snow and change temperature
        c[1] = cn * kb[2] / dzb[2]
        b[1] = c[1] + dzm[1] * crho[1] / dt
        d[1] = tg[1] * dzm[1] * crho[1] / dt +
            cn * gflux + 
            (1. - cn) * (kb[2] * (tg[2] - tg[1]) / dzb[2] + gflux0)
          
    # beginning of matrix calculation for a*T(z-1) - b*T(z) + c*T(z+1) + d  = 0

    tg[1] = b[1]

    for k = 2:km-1 # change to 2 to km if tgb is not used
        tg[k] = b[k] - a[k] / tg[k-1] * c[k-1]
        d[k] = d[k] + a[k] / tg[k-1] * d[k-1]
    end

    tg[km] = tgb

    for k = km-1:-1:1
        tg[k] = (d[k] + c[k] * tg[k+1]) / tg[k]
    end

    gflux = kb[1] * (tg0 - tg[1]) / dzb[1]

    return tg, gflux

end

function FRextended!(tgb::Float64, tgba::Float64, Dthermal::Float64, dzb::Float64, Tm::Float64)

    dt = 60. * 60. * 24. # 1 day
    omega = 2. * pi / (24. * 60. * 60. * 365.) # variation in the annual time scale
    Da = sqrt(2. * Dthermal / omega)
    C1 = 1. + 2. * dzb / Da

    tgb =
        (-2. * Dthermal / Da * (tgb - tgba) / dzb - omega * (tgb - Tm)) / C1 *
        dt + tgb


    return tgb
end

function soilinterp(zsoilobs::Vector{Float64}, tgini::Vector{Float64}, zbsoil::Vector{Float64},ngrid)
    # model grid
    zm = zeros(Float64, ngrid)
    zb = vcat(0., zbsoil)
    for k = 1:ngrid
        zm[k] = (zb[k] + zb[k+1]) * 0.5
    end
    itp = LinearInterpolation(zsoilobs, tgini, extrapolation_bc = Linear())
    tsoilini = itp(zm)
       return tsoilini   
end
