# flotech

`flotech` 
# Multiphase Flow Correlation Library

A Python library implementing various multiphase flow correlations for oil and gas production systems, including fluid property calculations and pressure gradient predictions.

## Features

- **Multiple Correlation Models**:
  - Mukherjee and Brill (1985)
  - Ansari et al. (1994)
  - Aziz and Govier (1972)
  - Stanford Correlation
- **Fluid Property Calculations**:
  - Oil, water, and gas properties
  - PVT (Pressure-Volume-Temperature) relationships
- **IPR Calculations**:
  - Vogel's Inflow Performance Relationship
  - Darcy's equation

## Installation


You can install `flotech` via pip:

```bash
pip install flotech

_______
## Usage


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


def pressure_traverse_general(oil_rate,corr_object):
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

def vlp(rates, corr_obj):
    bhps =[]
    for q in rates:
        p, dpdz = pressure_traverse_general(q,corr_obj )
        bhp = p[-1]
        bhps.append(bhp)
    return bhps

rates = np.linspace(10, 5000, sample_size)


bhps = vlp(rates,FB)
plt.plot(rates, bhps, label="FancherBrown"); 

bhps = vlp(rates,BB)
plt.plot(rates, bhps, label="Beggs and Brill"); 

bhps = vlp(rates,AZ)
plt.plot(rates, bhps, label="Aziz Govier"); 

bhps = vlp(rates,ST)
plt.plot(rates, bhps, label="Stanford University FM"); 

bhps = vlp(rates,HY)
plt.plot(rates, bhps, label="Hydro 3Phase"); 

bhps = vlp(rates,DR)
plt.plot(rates, bhps, label="Duns and Ross"); 

bhps = vlp(rates,MB)
plt.plot(rates, bhps, label="Mukherjee and Brill 1985"); 


bhps = vlp(rates,AN)
plt.plot(rates, bhps, label="Ansari 1994 - Mechanistic"); 


bhps = vlp(rates,GR)
plt.plot(rates, bhps, label="Gray 1974"); 

# Add this after all pressure traverse runs
# --- Vogel IPR parameters ---
pres = 5000  # psia (Reservoir pressure)
qmax = 5000  # STB/day (Maximum rate at zero pressure)

# Generate Vogel IPR data
pwfs = np.linspace(0, pres, 100)
q_vogel = qmax * (1 - 0.2 * (pwfs / pres) - 0.8 * (pwfs / pres) ** 2)


plt.plot(q_vogel, pwfs, label='Vogel IPR', color='black', linewidth=2)


pwfs = np.linspace(0, pres, 100)
q_vogel = 3500 * (1 - 0.2 * (pwfs / pres) - 0.8 * (pwfs / pres) ** 2)


plt.plot(q_vogel, pwfs, label='Vogel IPR', linewidth=2)


plt.xlabel('Oil Rate (STB/day)')
plt.ylabel('Bottomhole Pressure (psia)')
plt.title('Vogel IPR vs TPR')
plt.legend()
plt.grid(True)
plt.show()


plt.legend()
plt.show()

