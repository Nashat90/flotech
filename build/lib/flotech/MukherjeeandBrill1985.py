# coding=utf-8
import math
import flotech.FluidProps as FluidProps 

def Pgrad(P, T, oil_rate, wtr_rate, Gor, gas_grav, oil_grav, wtr_grav, d, angle):
    """
    Mukherjee and Brill (1985) pressure gradient correlation.
    Returns:
        dP_dz (float): Pressure gradient (psi/ft)
    """
    pi = math.pi
    g = 32.174  # ft/s²
    theta = angle * pi / 180
    Tsep = 60
    Psep = 114.7

    # --- Fluid properties ---
    Z = FluidProps.zfact((T + 460) / FluidProps.Tc(gas_grav), P / FluidProps.Pc(gas_grav))
    Wor = wtr_rate / oil_rate
    TDS = FluidProps.salinity(wtr_grav)
    Pb = FluidProps.Pbub(T, Tsep, Psep, gas_grav, oil_grav, Gor)
    Rso = FluidProps.sol_gor(T, P, Tsep, Psep, Pb, gas_grav, oil_grav)
    Rsw = FluidProps.sol_gwr(P, T, TDS)
    Bo = FluidProps.oil_fvf(T, P, Tsep, Psep, Pb, Rso, gas_grav, oil_grav)
    Bw = FluidProps.wtr_fvf(P, T, TDS)
    Bg = FluidProps.gas_fvf(P, T, gas_grav)
    muo = FluidProps.oil_visc(T, P, Tsep, Psep, Pb, Rso, gas_grav, oil_grav)
    muw = FluidProps.wtr_visc(P, T, TDS)
    mug = FluidProps.gvisc(P, T + 460, Z, gas_grav)
    rhoo = FluidProps.oil_dens(T, P, Tsep, Psep, Pb, Bo, Rso, gas_grav, oil_grav)
    rhow = 62.368 * wtr_grav / Bw
    rhog = 2.7 * gas_grav * P / (T + 460) / Z

    # --- Volumetric rates and velocities ---
    Axs = pi / 4 * (d / 12) ** 2
    qo = Bo * oil_rate / 15387
    qw = Bw * Wor * oil_rate / 15387
    ql = qo + qw
    qg = Bg * max(0, Gor - Rso - Rsw * Wor) * oil_rate / 86400
    vsl = ql / Axs
    vsg = qg / Axs
    vm = vsl + vsg

    # --- Flow regime (simplified approximation) ---
    if angle > 10:
        regime = "segregated"
    elif angle > 0:
        regime = "intermittent"
    else:
        regime = "distributed"

    # --- Liquid holdup correlation (Mukherjee & Brill) ---
    if regime == "segregated":
        hl = min(1.0, max(0.2, 1 - math.exp(-vsl / (vsg + 1e-6))))
    elif regime == "intermittent":
        hl = min(1.0, max(0.3, 0.8 * vsl / (vsg + vsl)))
    else:  # distributed
        hl = min(1.0, max(0.1, 0.5))

    # --- Mixture properties ---
    rhol = (Bw * Wor * rhow + Bo * rhoo) / (Bw * Wor + Bo)
    rhom = rhog * (1 - hl) + rhol * hl
    mul = (Bw * Wor * muw + Bo * muo) / (Bw * Wor + Bo)

    # --- Friction factor ---
    Nre = 1488 * rhom * vm * (d / 12) / mul
    fn = Fric(Nre, 0.0006)

    # --- Pressure gradient components ---
    dP_dz_ele = rhom * math.sin(theta) / 144
    dP_dz_fric = fn * rhom * vm**2 / (2 * (d / 12)) / 144
    dP_dz = dP_dz_ele + dP_dz_fric

    return dP_dz

def Fric(Nre, eps):
    """Chen’s friction factor (explicit)."""
    try:
        Temp = -4 * math.log10((eps / 3.7065) - (5.0452 / Nre) * 
               math.log10((eps**1.1098 / 2.8257) + (7.149 / Nre)**0.8981))
    except:
        return 0.005
    return (1 / Temp) ** 2

def Pwf_q_MB(FWHP, FWHT, Oil_Rate, Water_Rate, GOR, GasGrav, API, WaterGrav, ID, Angle, Depth, FBHT):
    """
    Pressure traverse using Mukherjee & Brill pressure gradient correlation.
    """
    DPs = []
    Temps = []
    PressList = []
    DepthList = []

    PressList.append(FWHP)
    Temps.append(FWHT)
    DepthList.append(0)
    DPs.append(0)

    nSteps = 60
    Tgrad = (FBHT - FWHT) / Depth

    for i in range(1, nSteps + 1):
        DeltaD = Depth / nSteps * i
        DepthList.append(DeltaD)
        T = Temps[i-1] + Tgrad * (DepthList[i] - DepthList[i-1])
        Temps.append(T)
        p = PressList[i-1] + DPs[i-1] * (DepthList[i] - DepthList[i-1])
        dp = Pgrad_MB(p, T, Oil_Rate, Water_Rate, GOR, GasGrav, API, WaterGrav, ID, Angle)
        DPs.append(dp)
        PressList.append(p)

    return PressList[-1]
