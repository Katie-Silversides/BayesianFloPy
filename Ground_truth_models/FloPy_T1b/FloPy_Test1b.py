# -*- coding: utf-8 -*-
"""
Created on Tue May 30 12:38:05 2023
First test of a flopy model for Modflow
@author: ksil8584
"""

from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import flopy
import os

## Set variables
# specified head
h1 = 100 

#Number of layers
Nlay = 10

#Number of rows and columns
N = 101

#lengths of the sides of the model
L = 400.0

#Aquifer thickness
H = 50.0

#Hydrolic conductivity
k = 1.0

#Well pumping rate
q = -1000

# recharge rate
rch_rate = 0

## Create FLoPy simulation object
# work_dir = 'C:\\Users\\ksil8584\\OneDrive - The University of Sydney (Staff)\\Documents\\Modflow\mf6.4.1\\examples\\FloPy_T1'
workspace = os.getcwd()
mf6_path = "C:\\Users\\ksil8584\\OneDrive - The University of Sydney (Staff)\\Documents\\Modflow\mf6.4.1\\bin\\mf6.exe"
name = "FloPy_Test1"
sim = flopy.mf6.MFSimulation(
    sim_name=name, exe_name=mf6_path, version="mf6", sim_ws=workspace
)

## Create FloPy temporal discretization (TDIS) object
tdis = flopy.mf6.ModflowTdis(
    sim, pname="tdis", time_units="DAYS", nper=1, perioddata=[(1.0, 1, 1.0)] # add  for _ in range(2)  to persiod data if nper >1
)

## Create FloPy iterative model solution (IMS) object
ims = flopy.mf6.ModflowIms(
    sim,
    pname="ims",
    complexity="SIMPLE",
    linear_acceleration="BICGSTAB",
)

## Create Flopy groundwater flow (gwf) model object
model_nam_file = f"{name}.nam"
gwf = flopy.mf6.ModflowGwf(
    sim,
    modelname=name,
    model_nam_file=model_nam_file,
    save_flows=True,
    newtonoptions="NEWTON UNDER_RELAXATION",
)

## Create discretization (DIS) Package
# All layers are given equal thickness
bot = np.linspace(-H / Nlay, -H, Nlay)
#  delrow and delcol are computed from model size L and number of cells
delrow = delcol = L / (N - 1)
dis = flopy.mf6.ModflowGwfdis(
    gwf,
    nlay=Nlay,
    nrow=N,
    ncol=N,
    delr=delrow,
    delc=delcol,
    top=0.0,
    botm=bot,
)

## Create initial conditions (IC) Package
start = h1 * np.ones((Nlay, N, N))
ic = flopy.mf6.ModflowGwfic(gwf, pname="ic", strt=start)

## Create node property flow (NPF) package
npf = flopy.mf6.ModflowGwfnpf(
    gwf,
    icelltype=1,
    k=k,
)
            
## Create storage (STO) package
sto = flopy.mf6.ModflowGwfsto(
    gwf, 
    iconvert = 1
)

## Create constant head (CHD) Package
chd_rec = []
layer = 0
for row_col in range(0, N):
    chd_rec.append(((layer, row_col, 0), h1))
    chd_rec.append(((layer, row_col, N - 1), h1))
    if row_col != 0 and row_col != N - 1:
        chd_rec.append(((layer, 0, row_col), h1))
        chd_rec.append(((layer, N - 1, row_col), h1))
chd = flopy.mf6.ModflowGwfchd(
    gwf,
    stress_period_data=chd_rec,
)

# To check the structured array for a particular persion (iper) use:
# iper = 0
# ra = chd.stress_period_data.get_data(key=iper)
# ra

## Create well (WEL) package
wel_rec = [(Nlay - 1, int(N / 4), int(N / 4), q)]
wel = flopy.mf6.ModflowGwfwel(
    gwf,
    stress_period_data=wel_rec,
)

#Creach recharge (rcha) package
rch = flopy.mf6.ModflowGwfrcha(
    gwf, 
    readasarrays=True, 
    pname="rch", 
    recharge=rch_rate
)

## Create output control (OC) Package
headfile = f"{name}.hds"
head_filerecord = [headfile]
budgetfile = f"{name}.cbb"
budget_filerecord = [budgetfile]
saverecord = [("HEAD", "ALL"), ("BUDGET", "ALL")]
printrecord = [("HEAD", "LAST")]
oc = flopy.mf6.ModflowGwfoc(
    gwf,
    saverecord=saverecord,
    head_filerecord=head_filerecord,
    budget_filerecord=budget_filerecord,
    printrecord=printrecord,
)

## Write the datasets
sim.write_simulation()

## Run simulation
success, buff = sim.run_simulation()
assert success, "MODFLOW 6 did not terminate normally."


### Post process head results
## Get data 
h = gwf.output.head().get_data(kstpkper=(0, 0))
x = y = np.linspace(0, L, N)
y = y[::-1]
vmin, vmax = 90.0, 110.0
contour_intervals = np.arange(90, 110.1, 1.0)

# create monitoring wells
m_well2 = np.empty([4,4])
for m_x in range(20,100,20):
    m_xi = int(m_x/20-1)
    for m_y in range(20,100,20):
        m_yi = int(m_y/20-1)
        m_well2[m_xi,m_yi] = h[0][m_x,m_y]
np.savetxt('FloPy_Test1b_monitoring_wells.csv',m_well2,delimiter=',')     

# create 100 monitoring wells 
m_well = np.empty([10,10])
for m_x in range(5,100,10):
    m_xi = int(m_x/10-0.5)
    for m_y in range(5,100,10):
        m_yi = int(m_y/10-0.5)
        m_well[m_xi,m_yi] = h[0][m_x,m_y]
# np.savetxt('FloPy_Test1b_monitoring_wells_100.csv',m_well,delimiter=',')        

# create 9x9 monitoring wells
m_well3 = np.empty([9,9])
for m_x in range(10,100,10):
    m_xi = int(m_x/10-1)
    for m_y in range(10,100,10):
        m_yi = int(m_y/10-1)
        m_well3[m_xi,m_yi] = h[0][m_x,m_y]
np.savetxt('FloPy_Test1b_monitoring_wells_9x9.csv',m_well3,delimiter=',')

##Plot layer 1
fig = plt.figure(figsize=(6, 6))
ax = fig.add_subplot(1, 1, 1, aspect="equal")
plt.rcParams['font.size'] = 14
c = ax.contour(x, y, h[0], contour_intervals, colors="black")
plt.clabel(c, fmt="%2.1f")
plt.xlabel("X (m)", size=16)
plt.ylabel("Y (m)", size=16)
plt.title("Pumping = 0, Recharge = 0.0075", size=16)
plt.axis('equal')
ax.set_aspect('auto')
ax.set_xbound(0,400)
plt.xticks(np.arange(0,401,100))
plt.yticks(np.arange(0,401,100))
