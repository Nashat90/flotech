# coding=utf-8
import math
import flotech.FluidProps as FluidProps   # Identical to original imports

# No flow regimes are returned from this model
def Pgrad(P, T, oil_rate, wtr_rate, Gor, gas_grav, oil_grav, wtr_grav, d, angle):
    """
    Calculate pressure gradient using Fancher & Brown (1963) correlation.
    Parameters identical to Beggs & Brill.
    Returns:
        dP_dz (float): Pressure gradient (psi/ft)
    """
    # --- Constants ---
    pi = math.pi
    g = 32.174  # ft/s²
    Psep = 114.7  # psia
    Tsep = 50     # °F

    # --- Convert angle to radians ---
    theta = angle * pi / 180

    # --- Fluid Properties (same as Beggs & Brill) ---
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
    rhog = 2.699 * gas_grav * P / (T + 460) / Z

    # --- Mixture Properties ---
    rhol = (Bw * Wor * rhow + Bo * rhoo) / (Bw * Wor + Bo)
    mul = (Bw * Wor * muw + Bo * muo) / (Bw * Wor + Bo)

    # --- Superficial Velocities (ft/s) ---
    Axs = pi / 4 * (d / 12) ** 2
    qo = Bo * oil_rate / 15387
    qw = Bw * Wor * oil_rate / 15387
    ql = qo + qw
    qg = Bg * max(0, Gor - Rso - Rsw * Wor) * oil_rate / 86400
    vsl = ql / Axs
    vsg = qg / Axs
    vm = vsl + vsg

    # --- Fancher & Brown Holdup Correlation ---
    # Empirical slip velocity (ft/s)
    Vs = 0.8 * (g * (d / 12) * (rhol - rhog) / rhol) ** 0.5  # Slip velocity
    yl = vsl / (vsl + vsg + Vs)  # Liquid holdup
    yl = max(vsl / vm, min(yl, 1.0))  # Bounded between no-slip and 1.0

    # --- Mixture Density & Viscosity ---
    rhom = rhol * yl + rhog * (1 - yl)
    mum = mul * yl + mug * (1 - yl)

    # --- Friction Factor (Chen equation, same as Beggs & Brill) ---
    Nre = 1488 * rhom * vm * (d / 12) / mum
    fn = Fric(Nre, 0.0006)

    # --- Pressure Gradient (psi/ft) ---
    dP_dz_ele = rhom * math.sin(theta) / 144  # Elevation
    dP_dz_fric = fn * rhom * vm ** 2 / (2 * (d / 12)) / 144  # Friction
    dP_dz = dP_dz_ele + dP_dz_fric  # Total (neglect acceleration)

    return dP_dz
def Fric(Nre, eps):
    """Chen's explicit friction factor (identical to Beggs & Brill)."""
    try:
        Temp = -4 * math.log10((eps / 3.7065) - (5.0452 / Nre) * 
               math.log10((eps ** 1.1098 / 2.8257) + (7.149 / Nre) ** 0.8981))
    except:
        return 0.005  # Fallback value
    return (1 / Temp) ** 2
def Pwf_q(FWHP, FWHT, Oil_Rate, Water_Rate, GOR, GasGrav, API, WaterGrav, ID, Angle, Depth, FBHT):
    """Identical to Beggs & Brill's version but uses Fancher & Brown's Pgrad()."""
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
        dp = Pgrad(p, T, Oil_Rate, Water_Rate, GOR, GasGrav, API, WaterGrav, ID, Angle)
        DPs.append(dp)
        PressList.append(p)

    return PressList[-1]