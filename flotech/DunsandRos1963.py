# coding=utf-8
import math
import flotech.FluidProps as FluidProps 

def Pgrad(P, T, oil_rate, wtr_rate, Gor, gas_grav, oil_grav, wtr_grav, d, angle):
    """
    Calculate pressure gradient using Duns and Ros (1963) correlation.
    Assumes vertical well flow.
    Returns:
        dP_dz (float): Pressure gradient (psi/ft)
    """
    pi = math.pi
    g = 32.174  # ft/s²
    Psep = 114.7  # psia
    Tsep = 60     # °F

    theta = angle * pi / 180

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
    rhog = 2.699 * gas_grav * P / (T + 460) / Z

    # --- Superficial velocities ---
    Axs = pi / 4 * (d / 12) ** 2
    qo = Bo * oil_rate / 15387
    qw = Bw * Wor * oil_rate / 15387
    ql = qo + qw
    qg = Bg * max(0, Gor - Rso - Rsw * Wor) * oil_rate / 86400
    vsl = ql / Axs
    vsg = qg / Axs
    vm = vsl + vsg

    # --- Liquid holdup (Duns & Ros empirical estimation) ---
    # Flow regime decision based on superficial velocities
    if vsg < 0.5:
        hl = 1.0  # bubble or slug flow, all liquid assumed to be held up
    elif vsg < 10:
        hl = 0.7  # transition/slug flow
    else:
        hl = 0.3  # mist flow, little liquid holdup

    # --- Mixture properties ---
    rhol = (Bw * Wor * rhow + Bo * rhoo) / (Bw * Wor + Bo)
    rhom = rhog * (1 - hl) + rhol * hl
    mul = (Bw * Wor * muw + Bo * muo) / (Bw * Wor + Bo)

    # --- Reynolds number & friction factor ---
    Nre = 1488 * rhom * vm * (d / 12) / mul
    fn = Fric(Nre, 0.0006)

    # --- Pressure gradient components ---
    dP_dz_ele = rhom * math.sin(theta) / 144
    dP_dz_fric = fn * rhom * vm**2 / (2 * (d / 12)) / 144
    dP_dz = dP_dz_ele + dP_dz_fric

    return dP_dz

def Fric(Nre, eps):
    """Chen's explicit friction factor."""
    try:
        Temp = -4 * math.log10((eps / 3.7065) - (5.0452 / Nre) * 
               math.log10((eps**1.1098 / 2.8257) + (7.149 / Nre)**0.8981))
    except:
        return 0.005
    return (1 / Temp) ** 2

def Pwf_q_DunsRos(FWHP, FWHT, Oil_Rate, Water_Rate, GOR, GasGrav, API, WaterGrav, ID, Angle, Depth, FBHT):
    """
    Pressure traverse using Duns and Ros pressure gradient correlation.
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
        dp = Pgrad(p, T, Oil_Rate, Water_Rate, GOR, GasGrav, API, WaterGrav, ID, Angle)
        DPs.append(dp)
        PressList.append(p)

    return PressList[-1]

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm

def ____(usg):
    """
    Duns & Ros empirical regime based on vsg.
    Returns:
        1 = Bubble/Slug
        2 = Slug/Transition
        3 = Mist
    """
    if usg < 0.5:
        return 1
    elif usg < 10:
        return 2
    elif usg>10:
        print("Mist found")
        return 3

def plot_dunsros_flow_map():
    usl_range = np.logspace(-2, 1, 300)     # 0.01 to 10 ft/s
    usg_range = np.logspace(-2, 2, 300)     # 0.01 to 100 ft/s

    USL, USG = np.meshgrid(usl_range, usg_range)

    usl_flat = USL.flatten()
    usg_flat = USG.flatten()
    regime_flat = np.zeros_like(usl_flat, dtype=int)

    # Assign flow regime per point
    for i in range(len(usl_flat)):
        regime_flat[i] = ____(usg_flat[i])

    # Define colors
    color_map = {1: 'blue', 2: 'green', 3: 'red'}
    labels = {1: 'Bubble/Slug', 2: 'Slug/Transition', 3: 'Mist'}
    colors = [color_map[r] for r in regime_flat]

    # Plot
    plt.figure(figsize=(10, 7))
    scatter = plt.scatter(usl_flat, usg_flat, c=colors, alpha=0.7, s=10, edgecolors='none')

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Liquid Superficial Velocity (ft/s, log) ------')
    plt.ylabel('Gas Superficial Velocity (ft/s, log)')
    plt.title('Duns & Ros Flow Regime Map (Colored Scatter, log-log)')
    plt.grid(True, which='both', ls='--', lw=0.5)

    # Manual legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=color_map[i], label=labels[i]) for i in sorted(labels)]
    plt.legend(handles=legend_elements, loc='upper right')

    plt.tight_layout()
    plt.show()
