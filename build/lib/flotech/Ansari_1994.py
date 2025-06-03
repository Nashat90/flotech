import flotech.FluidProps as FluidProps 
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm

def Pgrad(P, T, oil_rate, wtr_rate, Gor, gas_grav, oil_grav, wtr_grav, d, angle):
    
    """
    Calculate pressure gradient using Ansari et al. (1994) correlation.
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
    _vsg = vsg
    _vsl = vsl
    # --- Flow Regime Identification (Ansari's Approach) ---
    # (1) Bubble/Slug Transition
    N_fr = vm ** 2 / (g * (d / 12))
    if N_fr < 0.5:
        regime = "bubble"
    elif 0.5 <= N_fr < 1.5:
        regime = "slug"
    else:
        regime = "annular"
    _vsg = vsg
    _vsl = vsl
    _regime = regime
    # --- Holdup Calculation (Ansari Empirical) ---
    if regime == "bubble":
        yl = 1 - 0.25 * (vsg / vm)
    elif regime == "slug":
        yl = 0.9 - 0.2 * (vsg / vm)
    else:  # annular
        yl = 0.8 - 0.3 * (vsg / vm)
    yl = max(vsl / vm, min(yl, 1.0))  # Bounds: no-slip ≤ yl ≤ 1.0

    # --- Mixture Density & Viscosity ---
    rhom = rhol * yl + rhog * (1 - yl)
    mum = mul * yl + mug * (1 - yl)

    # --- Friction Factor (Reuse Beggs & Brill's `Fric()`) ---
    Nre = 1488 * rhom * vm * (d / 12) / mum
    fn = Fric(Nre, 0.0006)

    # --- Pressure Gradient (psi/ft) ---
    dP_dz_ele = rhom * math.sin(theta) / 144  # Elevation
    dP_dz_fric = fn * rhom * vm ** 2 / (2 * (d / 12)) / 144  # Friction
    dP_dz = dP_dz_ele + dP_dz_fric  # Total (neglect acceleration)

    return dP_dz
def Fric(Nre, eps):
    """Calculate Fanning Friction Factor using the Chen Equation """
    try:
        math.log
        Temp = -4 * math.log10((eps / 3.7065) - (5.0452 / Nre) * math.log10((eps ** 1.1098 / 2.8257) + (7.149 / Nre) ** 0.8981) ) 
    except Exception as inst:
         print(type(inst))    # the exception instance
         print(inst.args)     # arguments stored in .args
         print(inst)       
    
    return (1 / Temp) ** 2
def flow_regime_contour_map(d=2.4414):
    g = 32.174  # ft/s²

    vsl_range = np.linspace(0.01, 10, 200)  # liquid superficial velocity
    vsg_range = np.linspace(0.01, 30, 200)  # gas superficial velocity

    regime_map = np.empty((len(vsg_range), len(vsl_range)), dtype=int)

    # Map regimes to integers for plotting
    regime_dict = {"bubble": 0, "slug": 1, "annular": 2}

    for i, vsg in enumerate(vsg_range):
        for j, vsl in enumerate(vsl_range):
            vm = vsl + vsg
            N_fr = vm ** 2 / (g * (d / 12))

            if N_fr < 0.5:
                regime_map[i, j] = regime_dict["bubble"]
            elif 0.5 <= N_fr < 1.5:
                regime_map[i, j] = regime_dict["slug"]
            else:
                regime_map[i, j] = regime_dict["annular"]

    VSL, VSG = np.meshgrid(vsl_range, vsg_range)

    # Define colormap and boundaries for contour
    cmap = ListedColormap(['blue', 'green', 'red'])
    bounds = [-0.5, 0.5, 1.5, 2.5]  # boundaries between regime codes
    norm = BoundaryNorm(bounds, cmap.N)

    plt.figure(figsize=(10, 7))
    cs = plt.contourf(VSL, VSG, regime_map, levels=[-0.5,0.5,1.5,2.5], colors=['blue','green','red'], alpha=0.6)

    # Add contour lines
    contour_lines = plt.contour(VSL, VSG, regime_map, levels=[0,1,2], colors='k', linewidths=1)

    plt.clabel(contour_lines, fmt={0:'Bubble', 1:'Slug', 2:'Annular'}, inline=True, fontsize=10)

    plt.xlabel('Liquid Superficial Velocity, v_sl (ft/s)')
    plt.ylabel('Gas Superficial Velocity, v_sg (ft/s)')
    plt.title('Flow Regime Contour Map based on Ansari et al. (1994)')

    from matplotlib.patches import Patch
    handles = [Patch(color='blue', label='Bubble'),
               Patch(color='green', label='Slug'),
               Patch(color='red', label='Annular')]
    plt.legend(handles=handles, loc='upper right')
    plt.yscale("log"); plt.xscale("log")
    plt.grid(True)
    plt.show()
    return plt