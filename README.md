
# flotech

`flotech` 
# Multiphase Flow Correlation Library

A Python library implementing various multiphase flow correlations for oil and gas production systems, including fluid property calculations and pressure gradient predictions.

## Features

- **Multiple Correlation Models**:
  - Mukherjee and Brill (1985)
  - Ansari et al. (1994)
  - Aziz and Govier (1972)
  - Hydro Model
  - Gray
  - Fancher and Brown
  - Duns and Ross
  - Stanford Correlation
  - Plus Some more
- **Fluid Property Calculations**:
  - Oil, water, and gas properties
  - PVT (Pressure-Volume-Temperature) relationships
- **IPR Calculations**:
  - Vogel's Inflow Performance Relationship
  - Darcy's equation

## Installation
You can install `flotech` via pip:

## Installation 
```bash
pip install flotech
```
# Well Performance Analysis with FloTech

This project calculates Vertical Lift Performance (VLP) and compares different multiphase flow correlations along with Vogel Inflow Performance Relationship (IPR).

```python
import numpy as np
import pandas as pd
import math
import flotech
import matplotlib.pyplot as plt
import flotech.FluidProps
import flotech.BeggsandBrill as BB
import flotech.Gray1974 as GR
import flotech.Ansari_1994 as AN
from flotech import HYDRO3Phase as HY, AzizandGovier1972 as AZ, FancherandBrown1963 as FB, Stanford as ST
from flotech import DunsandRos1963 as DR, MukherjeeandBrill1985 as MB

# Input parameters
oil_rate = 900
water_rate = 500
gor = 250
gas_grav = 0.65
oil_grav = 40
wtr_grav = 1.07
diameter = 3
angle = 90.0
thp = 100
tht = 100.0
twf = 150.0
depth = 8000
sample_size = 150

# Temperature gradient calculation
def temp_gradient(t0, t1, depth):
    if depth == 0:
        return 0
    else:
        return abs(t0 - t1) / depth    

t_grad = temp_gradient(tht, twf, depth)
depths = np.linspace(0, depth, sample_size)
temps = tht + t_grad * depths

# Pressure traverse function using a given correlation
def pressure_traverse_general(oil_rate, corr_object):
    p = []
    dpdz = []
    for i in range(len(depths)):
        if i == 0:
            p.append(thp)
        else:
            dz = (depths[i] - depths[i-1])
            pressure = p[i-1] + dz * dpdz[i-1]
            p.append(pressure)
        dpdz_step = corr_object.Pgrad(p[i], temps[i], oil_rate, water_rate, gor, gas_grav, oil_grav, wtr_grav, diameter, angle) 
        dpdz.append(dpdz_step)
    return p, dpdz

# Calculate bottom-hole pressures (BHP) for a range of rates
def vlp(rates, corr_obj):
    bhps = []
    for q in rates:
        p, dpdz = pressure_traverse_general(q, corr_obj)
        bhp = p[-1]
        bhps.append(bhp)
    return bhps

rates = np.linspace(10, 5000, sample_size)

# Plot VLP curves for different correlations
bhps = vlp(rates, FB)
plt.plot(rates, bhps, label="FancherBrown")

bhps = vlp(rates, BB)
plt.plot(rates, bhps, label="Beggs and Brill")

bhps = vlp(rates, AZ)
plt.plot(rates, bhps, label="Aziz Govier")

bhps = vlp(rates, ST)
plt.plot(rates, bhps, label="Stanford University FM")

bhps = vlp(rates, HY)
plt.plot(rates, bhps, label="Hydro 3Phase")

bhps = vlp(rates, DR)
plt.plot(rates, bhps, label="Duns and Ross")

bhps = vlp(rates, MB)
plt.plot(rates, bhps, label="Mukherjee and Brill 1985")

bhps = vlp(rates, AN)
plt.plot(rates, bhps, label="Ansari 1994 - Mechanistic")

bhps = vlp(rates, GR)
plt.plot(rates, bhps, label="Gray 1974")

# Vogel IPR curve parameters and plot
pres = 5000  # psia (Reservoir pressure)
qmax = 5000  # STB/day (Maximum rate at zero pressure)

pwfs = np.linspace(0, pres, 100)
q_vogel = qmax * (1 - 0.2 * (pwfs / pres) - 0.8 * (pwfs / pres) ** 2)

plt.plot(q_vogel, pwfs, label='Vogel IPR', color='black', linewidth=2)

plt.xlabel('Oil Rate (STB/day)')
plt.ylabel('Bottomhole Pressure (psia)')
plt.title('Vogel IPR vs TPR')
plt.legend()
plt.grid(True)
plt.show()

```


![Well Diagram](images/Full_System.png)


# Wellbore Pressure Profile 

This project calculates Vertical Lift Performance (VLP) and compares different multiphase flow correlations along the wellbore path


```python
import numpy as np
import pandas as pd
import math
import flotech
import matplotlib.pyplot as plt
import flotech.FluidProps
import flotech.BeggsandBrill as BB
import flotech.Gray1974 as GR
import flotech.Ansari_1994 as AN
from flotech import HYDRO3Phase as HY, AzizandGovier1972 as AZ, FancherandBrown1963 as FB, Stanford as ST
from flotech import DunsandRos1963 as DR, MukherjeeandBrill1985 as MB
import numpy as np

plt.figure(figsize=(10,20))
oil_rate= 900
water_rate = 500
gor = 250
gas_grav = 0.65
oil_grav = 40
wtr_grav = 1.07
diameter=3
angle=90.0
thp= 100
tht=100.0
twf=150.0
depth = 8000
sample_size =150

def temp_gradient(t0,t1, depth):
    if depth==0:
        return 0
    else:
        return abs(t0-t1)/depth    

t_grad = temp_gradient(tht,twf, depth)
depths = np.linspace(0, depth, sample_size)
temps = tht + t_grad * depths


def pressure_traverse_general(corr_object):
    p=[]
    dpdz=[]
    for i in range(len(depths)):

        if i==0:
            p.append(thp)
        else:
            dz = (depths[i]-depths[i-1])
            pressure = p[i-1]+dz*dpdz[i-1]
            p.append(pressure)

        dpdz_step = corr_object.Pgrad(p[i], temps[i], oil_rate, water_rate, gor, gas_grav, oil_grav, wtr_grav, diameter, angle) 
        dpdz.append(dpdz_step)
    return p, dpdz


bhps,_ = pressure_traverse_general(FB)
plt.plot(bhps, depths, label="FancherBrown", lw=3); 

bhps,_ = pressure_traverse_general(BB)
plt.plot(bhps, depths, label="Beggs and Brill", lw=3); 

bhps,_ = pressure_traverse_general(AZ)
plt.plot(bhps, depths, label="Aziz Govier", lw=3); 

bhps,_ = pressure_traverse_general(ST)
plt.plot(bhps, depths, label="Stanford University FM", lw=3); 

bhps,_ = pressure_traverse_general(HY)
plt.plot(bhps, depths, label="Hydro 3Phase", lw=3); 

bhps,_ = pressure_traverse_general(DR)
plt.plot(bhps, depths, label="Duns and Ross", lw=3); 

bhps,_ = pressure_traverse_general(MB)
plt.plot(bhps, depths, label="Mukherjee and Brill 1985", lw=3); 


bhps,_ = pressure_traverse_general(AN)
plt.plot(bhps, depths, label="Ansari 1994 - Mechanistic", lw=3); 


bhps,_ = pressure_traverse_general(GR)
plt.plot(bhps, depths, label="Gray 1974", lw=3); 

plt.gca().invert_yaxis()
plt.ylabel("TVD ft")
plt.xlabel("Pressure psig")

# Move x-axis ticks and labels to the top
plt.gca().xaxis.set_ticks_position('top')      # Move ticks to top
plt.gca().xaxis.set_label_position('top')      # Move x-axis label to top
plt.gca().tick_params(axis='x', which='both', direction='inout')  # Optional: change tick direction
plt.grid(alpha=.2)


plt.legend()
plt.show()
plt.savefig(r"images\well_pressures.png")
```
---
# Wellbore Plot
![Well Diagram](images/well_pressures.png)
