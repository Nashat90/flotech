import math
import flotech.FluidProps as FluidProps 

def Pgrad(P, T, oil_rate, wtr_rate, Gor, gas_grav, oil_grav, wtr_grav, d, angle):
    """
    Function to Calculate the Flowing Pressure Gradient by the Method of Gray (1974)
    
    Parameters:
    P          pressure, psia
    T          temperature, °F
    oil_rate   oil flowrate, stb/d
    wtr_rate   water flowrate, stb/d
    Gor        producing gas-oil ratio, scf/stb
    gas_grav   gas specific gravity
    oil_grav   API oil gravity
    wtr_grav   water specific gravity
    d          pipe I.D., in.
    angle      angle of pipe inclination in degrees (90° = vertical)
    """
    
    # Constants
    pi = math.pi
    Psep = 114.7  # Separator pressure, psia
    Tsep = 50     # Separator temperature, °F
    g = 32.174    # gravitational acceleration, ft/s²
    
    # Convert angle to radians
    theta = angle * pi / 180
    
    # Calculate fluid properties (same as Beggs-Brill)
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
    
    # Volume fraction weighted liquid properties
    rhol = (Bw * Wor * rhow + Bo * rhoo) / (Bw * Wor + Bo)
    mul = (Bw * Wor * muw + Bo * muo) / (Bw * Wor + Bo)
    
    # Calculate flow rates (ft³/s)
    qo = Bo * oil_rate / 15387
    qw = Bw * Wor * oil_rate / 15387
    ql = qo + qw
    qg = Bg * max(0, (Gor - Rso - Rsw * Wor)) * oil_rate / 86400
    
    # Superficial velocities (ft/s)
    Axs = pi / 4 * (d / 12) ** 2
    vsl = ql / Axs
    vsg = qg / Axs
    vm = vsl + vsg
    
    # Gas-liquid ratio
    GL_R = qg / ql if ql > 0 else 0
    
    # Calculate liquid holdup using Gray correlation
    # Gray's empirical correlation for holdup
    C1 = 1.0714 - (0.2218 * (vm)**2 / (d/12))
    C1 = max(0.13, min(C1, 1.0))
    
    yl = C1 - (0.0212 + 0.00423 * (vm)) * (GL_R)**0.5
    yl = max(vsl/vm, min(yl, 1.0))  # Ensure holdup is between no-slip and 1
    
    # Mixture density
    rhom = rhol * yl + rhog * (1 - yl)
    
    # Mixture viscosity for friction calculation
    mum = mul * yl + mug * (1 - yl)
    
    # Reynolds number
    Nre = 1488 * rhom * vm * (d / 12) / mum
    
    # Friction factor (using Chen equation as in Beggs-Brill)
    fn = Fric(Nre, 0.0006)
    
    # Pressure gradient components
    # Elevation component
    dP_dz_ele = rhom * math.sin(theta) / 144
    
    # Friction component
    dP_dz_fric = fn * rhom * vm**2 / (2 * (d/12)) / 144
    
    # Acceleration component (usually small in vertical flow)
    dP_dz_acc = 0  # Often neglected in Gray's method
    
    # Total pressure gradient (psi/ft)
    dP_dz = dP_dz_ele + dP_dz_fric + dP_dz_acc
    
    return dP_dz

# Reuse the same Fric function from Beggs and Brill implementation
def Fric(Nre, eps):
    """Calculate Fanning Friction Factor using the Chen Equation"""
    try:
        Temp = -4 * math.log10((eps / 3.7065) - (5.0452 / Nre) * 
               math.log10((eps ** 1.1098 / 2.8257) + (7.149 / Nre) ** 0.8981))
    except Exception as inst:
        print(type(inst))    # the exception instance
        print(inst.args)     # arguments stored in .args
        print(inst)       
        return 0.005  # default friction factor if calculation fails
    
    return (1 / Temp) ** 2

def Pwf_q_Gray(FWHP, FWHT, Oil_Rate, Water_Rate, GOR, GasGrav, API, WaterGrav, ID, Angle, Depth, FBHT):
    """
    Function to calculate the bottomhole flowing pressure (Pwf) using Gray correlation
    Similar to the Beggs-Brill version but using Gray's method for pressure gradient
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
        dp = Pgrad_Gray(p, T, Oil_Rate, Water_Rate, GOR, GasGrav, API, WaterGrav, ID, Angle)
        DPs.append(dp)
        
        PressList.append(p)
    
    return PressList[-1]