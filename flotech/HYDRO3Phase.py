# coding=utf-8
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
import math
import flotech.FluidProps as FluidProps  # Assuming same fluid properties module
 
def Pgrad(P, T, oil_rate, wtr_rate, Gor, gas_grav, oil_grav, wtr_grav, d, angle):
    """
    Calculate pressure gradient using the HYDRO correlation for multiphase flow.
    
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

    # Convert angle to radians
    theta = math.radians(angle)

    # Fluid properties (using same helper functions as Beggs & Brill)
    Z = FluidProps.zfact((T + 460) / FluidProps.Tc(gas_grav), P / FluidProps.Pc(gas_grav))
    Wor = wtr_rate / oil_rate
    TDS = FluidProps.salinity(wtr_grav)
    Pb = FluidProps.Pbub(T, Tsep=50, Psep=114.7, gas_grav=gas_grav, oil_grav=oil_grav, Gor=Gor)
    Rso = FluidProps.sol_gor(T, P, Tsep=50, Psep=114.7, Pb=Pb, gas_grav=gas_grav, oil_grav=oil_grav)
    Bo = FluidProps.oil_fvf(T, P, Tsep=50, Psep=114.7, Pb=Pb,Rs=Rso, gas_grav=gas_grav, oil_grav=oil_grav)
    Bw = FluidProps.wtr_fvf(P, T, TDS)
    Bg = FluidProps.gas_fvf(P, T, gas_grav)
    rho_o = FluidProps.oil_dens(T, P, Tsep=50, Psep=114.7, Pb=Pb, Bo=Bo, Rs=Rso, gas_grav=gas_grav, oil_grav=oil_grav)
    rho_w = 62.368 * wtr_grav / Bw
    rho_g = 2.699 * gas_grav * P / (T + 460) / Z

    # Mixture properties
    rho_l = (Bo * rho_o + Bw * Wor * rho_w) / (Bo + Bw * Wor)  # Liquid density (lb/ft³)
    q_l = (Bo * oil_rate + Bw * wtr_rate) / 15387  # Liquid flow rate (ft³/s)
    q_g = Bg * (Gor - Rso) * oil_rate / 86400 if (Gor - Rso) > 0 else 0  # Gas flow rate (ft³/s)
    A = pi * (d / 12)**2 / 4  # Pipe cross-sectional area (ft²)

    # Superficial velocities
    u_sl = q_l / A  # Liquid superficial velocity (ft/s)
    u_sg = q_g / A  # Gas superficial velocity (ft/s)
    u_m = u_sl + u_sg  # Mixture velocity (ft/s)

    # Flow regime determination (HYDRO-specific)
    regime = HYDRO_flow_regime(u_sl, u_sg, d, rho_l, rho_g)

    # Liquid holdup (HYDRO model)
    H_l = HYDRO_holdup(regime, u_sl, u_sg, d, theta, rho_l, rho_g)

    # Mixture density
    rho_m = H_l * rho_l + (1 - H_l) * rho_g

    # Friction factor (Chen equation, same as Beggs & Brill)
    mu_l = FluidProps.oil_visc(T, P, Tsep=50, Psep=114.7, Pb=Pb, Rs=Rso, gas_grav=gas_grav, oil_grav=oil_grav)
    mu_g = FluidProps.gvisc(P, T + 460, Z, gas_grav)
    mu_m = H_l * mu_l + (1 - H_l) * mu_g
    Re = 1488 * rho_m * u_m * (d / 12) / mu_m
    f = Fric(Re, 0.0006)  # Using same friction factor function

    # Pressure gradient components
    dP_dz_grav =  0.052*(rho_m/7.4805) # Gravitational (psi/ft)
    dP_dz_fric = 2 * f * rho_m * u_m**2 / (d / 12) / 144  # Frictional (psi/ft)

    return dP_dz_grav + dP_dz_fric

def HYDRO_flow_regime(u_sl, u_sg, d, rho_l, rho_g):
    """
    Determine flow regime for HYDRO correlation.
    Returns:
        1 = Bubble
        2 = Slug
        3 = Stratified
        4 = Annular
        5 = Mist
    """
    u_m = u_sl + u_sg
    lambda_l = u_sl / u_m  # No-slip liquid holdup

    # HYDRO-specific regime boundaries
    if u_sg < 0.3 and lambda_l > 0.7:
        return 1  # Bubble
    elif u_sg < 10 and lambda_l > 0.3:
        return 2  # Slug
    elif u_m < 5 and lambda_l < 0.5:
        return 3  # Stratified
    elif u_sg > 20 and lambda_l < 0.2:
        return 4  # Annular
    else:
        return 5  # Mist

def HYDRO_holdup(regime, u_sl, u_sg, d, theta, rho_l, rho_g):
    """
    Calculate liquid holdup for HYDRO correlation.
    """
    u_m = u_sl + u_sg
    lambda_l = u_sl / u_m

    # HYDRO holdup model coefficients (empirical)
    if regime == 1:  # Bubble
        C = 0.85
    elif regime == 2:  # Slug
        C = 0.75
    elif regime == 3:  # Stratified
        C = 0.5
    elif regime == 4:  # Annular
        C = 0.1
    else:  # Mist
        C = 0.05

    # Inclination correction
    H_l = C * lambda_l * (1 + 0.1 * math.sin(theta))
    return min(max(H_l, lambda_l), 1.0)  # Clamped to [lambda_l, 1.0]

# Reuse Fric() from Beggs & Brill
def Fric(Nre, eps):
    
    """Chen's friction factor equation."""
    try:
        temp = -4 * math.log10(eps / 3.7065 - 5.0452 / Nre * math.log10(
            eps**1.1098 / 2.8257 + (7.149 / Nre)**0.8981))
        return (1 / temp)**2
    except:
        return 0.005  # Default if calculation fails

def plot_hydro_flow_map(d=2.4414, rho_l=50, rho_g=0.1):
    # Define velocity ranges
    u_sl_range = np.linspace(0.01, 10, 200)
    u_sg_range = np.linspace(0.01, 30, 200)
    USL, USG = np.meshgrid(u_sl_range, u_sg_range)

    regime_map = np.zeros_like(USL, dtype=int)

    # Populate regime map
    for i in range(USL.shape[0]):
        for j in range(USL.shape[1]):
            regime_map[i, j] = HYDRO_flow_regime(USL[i, j], USG[i, j], d, rho_l, rho_g)

    # Plotting
    cmap = ListedColormap(['orange', 'blue', 'green', 'orange', 'red', 'purple'])
    bounds = [0.5, 1.5, 2.5, 3.5, 4.5, 5.5]
    norm = BoundaryNorm(bounds, cmap.N)

    plt.figure(figsize=(10, 7))
    contour = plt.contourf(USL, USG, regime_map, levels=[1,2,3,4,5,6], 
                           cmap=cmap, norm=norm, alpha=.7)
    cbar = plt.colorbar(contour, ticks=[1,2,3,4,5])
    cbar.ax.set_yticklabels(['Bubble', 'Slug', 'Stratified', 'Annular', 'Mist'])

    plt.xlabel('Liquid Superficial Velocity (ft/s)')
    plt.ylabel('Gas Superficial Velocity (ft/s)')
    plt.title('HYDRO Flow Regime Contour Map')
    plt.xscale("log");plt.yscale("log");
    plt.grid(True)
    plt.show()
    