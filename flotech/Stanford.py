# coding=utf-8
import math
import flotech.FluidProps as FluidProps   # Reusing the same fluid properties module

def Pgrad(P, T, oil_rate, wtr_rate, Gor, gas_grav, oil_grav, wtr_grav, d, angle):
    """
    Calculate pressure gradient using the Stanford Flow Correlation.
    
    Parameters:
        P (float): Pressure (psia)
        T (float): Temperature (°F)
        oil_rate (float): Oil flow rate (stb/d)
        wtr_rate (float): Water flow rate (stb/d)
        Gor (float): Gas-oil ratio (scf/stb)
        gas_grav (float): Gas specific gravity (air=1.0)
        oil_grav (float): API oil gravity
        wtr_grav (float): Water specific gravity
        d (float): Pipe inner diameter (inches)
        angle (float): Pipe inclination angle (degrees, 0°=horizontal, 90°=vertical)
    
    Returns:
        float: Pressure gradient (psi/ft)
    """
    # Constants
    g = 32.174  # Gravitational acceleration (ft/s²)
    pi = math.pi
    g = 32.174  # ft/s²
    Psep = 114.7  # psia
    Tsep = 50     # °F
    # Convert angle to radians
    theta = math.radians(angle)

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
    rho_o = FluidProps.oil_dens(T, P, Tsep, Psep, Pb, Bo, Rso, gas_grav, oil_grav)
    rho_w = 62.368 * wtr_grav / Bw
    rho_g = 2.699 * gas_grav * P / (T + 460) / Z

    # --- Mixture Properties ---
    rhol = (Bw * Wor * rho_w + Bo * rho_o) / (Bw * Wor + Bo)
    mul = (Bw * Wor * muw + Bo * muo) / (Bw * Wor + Bo)

    # Mixture properties
    rho_l = (Bo * rho_o + Bw * Wor * rho_w) / (Bo + Bw * Wor)  # Liquid density (lb/ft³)
    q_l = (Bo * oil_rate + Bw * wtr_rate) / 15387  # Liquid flow rate (ft³/s)
    q_g = Bg * (Gor - Rso) * oil_rate / 86400 if (Gor - Rso) > 0 else 0  # Gas flow rate (ft³/s)
    A = pi * (d / 12)**2 / 4  # Pipe cross-sectional area (ft²)

    # Superficial velocities
    u_sl = q_l / A  # Liquid (ft/s)
    u_sg = q_g / A  # Gas (ft/s)
    u_m = u_sl + u_sg  # Mixture (ft/s)

    # Flow regime determination (Stanford-specific)
    regime = Stanford_flow_regime(u_sl, u_sg, d, rho_l, rho_g)

    # Liquid holdup (Stanford model)
    H_l = Stanford_holdup(regime, u_sl, u_sg, d, theta, rho_l, rho_g)

    # Mixture density
    rho_m = H_l * rho_l + (1 - H_l) * rho_g

    # Friction factor (Stanford uses a modified Colebrook-White)
    mu_l = FluidProps.oil_visc(T, P, Tsep=50, Psep=114.7, Pb=Pb, Rs=Rso, gas_grav=gas_grav, oil_grav=oil_grav)
    mu_g = FluidProps.gvisc(P, T + 460, Z, gas_grav)
    mu_m = H_l * mu_l + (1 - H_l) * mu_g
    Re = 1488 * rho_m * u_m * (d / 12) / mu_m
    f = Stanford_fric(Re, d)

    # Pressure gradient components
    dP_dz_grav = 0.052*(rho_m/7.4805) # Gravitational (psi/ft)
    dP_dz_fric = f * rho_m * u_m**2 / (2 * (d / 12)) / 144  # Frictional (psi/ft)

    return dP_dz_grav + dP_dz_fric

def Stanford_flow_regime(u_sl, u_sg, d, rho_l, rho_g):
    """
    Stanford flow regime classification.
    Returns:
        1 = Bubble
        2 = Slug
        3 = Churn
        4 = Annular
    """
    u_m = u_sl + u_sg
    lambda_l = u_sl / u_m  # No-slip liquid holdup

    # Stanford regime boundaries (empirical)
    if u_sg < 0.5 and lambda_l > 0.8:
        return 1  # Bubble
    elif u_sg < 10 and lambda_l > 0.25:
        return 2  # Slug
    elif u_m < 20 and lambda_l < 0.7:
        return 3  # Churn
    else:
        return 4  # Annular

def Stanford_holdup(regime, u_sl, u_sg, d, theta, rho_l, rho_g):
    """
    Stanford liquid holdup model.
    """
    u_m = u_sl + u_sg
    lambda_l = u_sl / u_m

    # Stanford coefficients
    if regime == 1:  # Bubble
        H_l = 0.98 * lambda_l**0.7
    elif regime == 2:  # Slug
        H_l = 0.83 * lambda_l**0.5
    elif regime == 3:  # Churn
        H_l = 0.65 * lambda_l**0.3
    else:  # Annular
        H_l = 0.1 * lambda_l**0.1

    # Inclination correction
    H_l *= (1 + 0.2 * math.sin(theta))
    return min(max(H_l, lambda_l), 1.0)  # Clamped to [lambda_l, 1.0]

def Stanford_fric(Re, d):
    """
    Stanford-modified Colebrook-White friction factor.
    """
    roughness = 0.0006  # Pipe roughness (ft)
    try:
        f = 1.02 / (4 * math.log10(roughness / (3.7 * (d / 12)) + 5.3 / Re**0.9))**2
        return f
    except:
        return 0.005  # Fallback value