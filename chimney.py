print("============================================")
print("Start Chimney EQ ground motion with gravity ")

import openseespy.postprocessing.ops_vis as opsv
import openseespy.postprocessing.Get_Rendering as opsplt
from openseespy.opensees import *
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
import ReadRecord
from PEERGM import processNGAfile

# SET UP ---------------------------------------------------------------------------------------------
wipe()  # clear opensees model
model('basic', '-ndm', 2, '-ndf', 3)  # 2 dimensions, 3 dof per node

# ----------------------------------------------------------------------------------------------------
# ----------------------- CONCRETE C40 ---------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------
IDconcCoreC40 = 1  # material ID tag -- confined core concrete
IDconcCoverC40 = 2  # material ID tag -- unconfined cover concrete
# nominal concrete compressive strength
fcC40 = -26.8  # CONCRETE Compressive Strength, MPa   (+Tension, -Compression)
EcC40 = 32500  # Concrete Elastic Modulus
# confined concrete
KfcC40 = 1.3  # ratio of confined to unconfined concrete strength
fc1CC40 = KfcC40 * fcC40  # CONFINED concrete (mander model), maximum stress
eps1CC40 = 2 * fc1CC40 / EcC40  # strain at maximum stress
fc2CC40 = 0.2 * fc1CC40  # ultimate stress
eps2CC40 = 5 * eps1CC40  # strain at ultimate stress
# unconfined concrete
fc1UC40 = fcC40  # UNCONFINED concrete (todeschini parabolic model), maximum stress
eps1UC40 = -0.003  # strain at maximum strength of unconfined concrete
fc2UC40 = 0.2 * fc1UC40  # ultimate stress
eps2UC40 = -0.01  # strain at ultimate stress
lamC40 = 0.1  # ratio between unloading slope at $eps2 and initial slope $Ec
# tensile-strength properties
ftCC40 = -0.14 * fc1CC40  # tensile strength +tension
ftUC40 = -0.14 * fc1UC40  # tensile strength +tension
EtsC40 = ftUC40 / 0.002  # tension softening stiffness

# ----------------------------------------------------------------------------------------------------
# ----------------------- CONCRETE C30 ---------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------
IDconcCoreC30 = 3  # material ID tag -- confined core concrete
IDconcCoverC30 = 4  # material ID tag -- unconfined cover concrete
# nominal concrete compressive strength
fcC30 = -14.3  # CONCRETE Compressive Strength, MPa   (+Tension, -Compression)
EcC30 = 30000  # Concrete Elastic Modulus
# confined concrete
KfcC30 = 1.3  # ratio of confined to unconfined concrete strength
fc1CC30 = KfcC30 * fcC30  # CONFINED concrete (mander model), maximum stress
eps1CC30 = 2 * fc1CC30 / EcC30  # strain at maximum stress
fc2CC30 = 0.2 * fc1CC30  # ultimate stress
eps2CC30 = 5 * eps1CC30  # strain at ultimate stress
# unconfined concrete
fc1UC30 = fcC30  # UNCONFINED concrete (todeschini parabolic model), maximum stress
eps1UC30 = -0.003  # strain at maximum strength of unconfined concrete
fc2UC30 = 0.2 * fc1UC30  # ultimate stress
eps2UC30 = -0.01  # strain at ultimate stress
lamC30 = 0.1  # ratio between unloading slope at $eps2 and initial slope $Ec
# tensile-strength properties
ftCC30 = -0.14 * fc1CC30  # tensile strength +tension
ftUC30 = -0.14 * fc1UC30  # tensile strength +tension
EtsC30 = ftUC30 / 0.002  # tension softening stiffness

# ----------------------------------------------------------------------------------------------------
# ----------------------- STEEL ----------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------
IDreinf = 5  # material ID tag -- reinforcement
Fy = 360  # STEEL yield stress
Es = 200000  # modulus of steel
Bs = 0.01  # strain-hardening ratio
R0 = 18.5  # control the transition from elastic to plastic branches
cR1 = 0.925  # control the transition from elastic to plastic branches
cR2 = 0.15  # control the transition from elastic to plastic branches

# Define materials for nonlinear columns
# -------------------------------------------------------------------------------------------------
# CORE CONCRETE C40 (confined)
uniaxialMaterial('Concrete02', IDconcCoreC40, fc1CC40, eps1CC40, fc2CC40, eps2CC40, lamC40, ftCC40, EtsC40)
# COVER CONCRETE C40 (unconfined)
uniaxialMaterial('Concrete02', IDconcCoverC40, fc1UC40, eps1UC40, fc2UC40, eps2UC40, lamC40, ftUC40, EtsC40)
# ------------------------------------------------------------------------------------------------
# CORE CONCRETE C30 (confined)
uniaxialMaterial('Concrete02', IDconcCoreC30, fc1CC30, eps1CC30, fc2CC30, eps2CC30, lamC30, ftCC30, EtsC30)
# COVER CONCRETE C30 (unconfined)
uniaxialMaterial('Concrete02', IDconcCoverC30, fc1UC30, eps1UC30, fc2UC30, eps2UC30, lamC30, ftUC30, EtsC30)
# ------------------------------------------------------------------------------------------------
# STEEL Reinforcing steel
uniaxialMaterial('Steel02', IDreinf, Fy, Es, Bs, R0, cR1, cR2)
# END MATERIAL parameters ------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------
# ----------------------- Different GEOMETRY Coordinate------------- Unit : m ------------------------------------
# ----------------------------------------------------------------------------------------------------------------
Unite = 10 ** 3  # ----- m to mm
#                             1   2   3   4   5   6   7   8   9   10  11  12  13  14  15  16  17  18  19  20  21  22  23  24
element_length = np.array(
    [0, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
     4]) * Unite  # Number of element 25
# print(element_length)
#                              1  2   3   4   5   6   7   8   9   10  11   12   13   14   15   16   17   18   19   20   21   22   23   24   25
Nodes_coordinate_Y = np.array(
    [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230,
     234]) * Unite
# Section Diameter ---------------         1     2     3     4     5     6     7     8     9     10    11    12    13    14    15    16    17    18    19    20    21    22    23    24
DSec = Column_Section_Diameter = np.array(
    [0, 27.312, 26.712, 26.112, 25.512, 24.912, 24.414, 23.740, 23.556, 23.156, 22.756, 22.356, 22.156, 21.956, 21.756,
     21.556, 21.356, 21.156, 20.956, 20.956, 20.956, 20.956, 20.956, 20.956, 20.956]) * Unite
# Section Thickness ---------------         1     2     3     4     5     6     7     8     9     10    11    12    13    14    15    16    17    18    19    20    21    22    23    24
TSec = Column_Section_Thickness = np.array(
    [0, 0.75, 0.73, 0.71, 0.69, 0.67, 0.65, 0.65, 0.65, 0.63, 0.60, 0.57, 0.55, 0.53, 0.51, 0.49, 0.45, 0.43, 0.41,
     0.39, 0.39, 0.39, 0.35, 0.35, 0.35]) * Unite
# --------------------------------------------
# Vertical reinforcement on the outer side of the tube wall (HRB400 grade is adopted for the reinforcement)
Rebar_spacing_VO = np.array(
    [0, 0.205, 0.205, 0.205, 0.215, 0.215, 0.225, 0.225, 0.225, 0.225, 0.235, 0.235, 0.235, 0.235, 0.235, 0.235, 0.235,
     0.235, 0.235, 0.235, 0.235, 0.235, 0.235, 0.235, 0.235, 0.235, 0.235, 0.235, 0.235]) * Unite
Rebar_diameter_VO = np.array(
    [0, 0.022, 0.022, 0.022, 0.022, 0.022, 0.022, 0.022, 0.022, 0.022, 0.022, 0.022, 0.022, 0.022, 0.022, 0.022, 0.022,
     0.020, 0.020, 0.020, 0.020, 0.020, 0.018, 0.018, 0.018, 0.018, 0.018, 0.018, 0.018]) * Unite
# Vertical reinforcement on the inner side of the tube wall (HRB400 grade is adopted for the reinforcement)
Rebar_spacing_VI = np.full(len(element_length), 0.200) * Unite
Rebar_diameter_VI = np.array(
    [0, 0.022, 0.022, 0.022, 0.022, 0.022, 0.022, 0.022, 0.022, 0.022, 0.022, 0.022, 0.022, 0.022, 0.022, 0.022, 0.022,
     0.020, 0.020, 0.020, 0.020, 0.020, 0.018, 0.018, 0.018, 0.018, 0.018, 0.018, 0.018]) * Unite

# Circumferential reinforcement on the outside of the wall (HRB400 grade is adopted for the reinforcement)
# Rebar_spacing_CO = np.full(len(element_length), 0.200)
# Rebar_diameter_CO = np.array([0.022, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.018, 0.018, 0.018, 0.018, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016])
# Circumferential steel bars on the inner side of the tube wall (HRB400 grade is adopted for the reinforcement)
# Rebar_spacing_CI = np.array([0.200, 0.200, 0.200, 0.200, 0.200, 0.200, 0.200, 0.200, 0.200, 0.200, 0.200, 0.200, 0.200, 0.200, 0.200, 0.200, 0.200, 0.200, 0.200, 0.200, 0.200, 0.200, 0.190, 0.175, 0.165, 0.140, 0.135, 0.135])
# Rebar_diameter_CI = np.array([0.018, 0.018, 0.018, 0.018, 0.018, 0.018, 0.018, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014])

Circumference_Section_Out = np.zeros(len(element_length))
Circumference_Section_In = np.zeros(len(element_length))
numBarsSec_Out = Number_Reinforcement_Out = np.zeros(len(element_length), dtype='i')
numBarsSec_In = Number_Reinforcement_In = np.zeros(len(element_length), dtype='i')
barAreaSec_Out = Area_Reinforcement_Out = np.zeros(len(element_length))
barAreaSec_In = Area_Reinforcement_In = np.zeros(len(element_length))
ColSecTag = Column_Section_Tag = np.zeros(len(element_length), dtype='i')
ro = np.zeros(len(element_length))
rs = np.zeros(len(element_length))
surf = np.zeros(len(element_length))
ri = np.zeros(len(element_length))
rc = np.zeros(len(element_length))
rc2 = np.zeros(len(element_length))
theta_Out = np.zeros(len(element_length))
theta_In = np.zeros(len(element_length))
nfCoreT = np.zeros(len(element_length), dtype='i')
nfCoverT = np.zeros(len(element_length), dtype='i')

CoverSec = 0.05 * Unite
for i in range(1, len(element_length)):
    Circumference_Section_Out[i] = (DSec[i] - (2 * CoverSec)) * math.pi  # Unit : m,
    Circumference_Section_In[i] = ((DSec[i] - (2 * TSec[i])) + (2 * CoverSec)) * math.pi  # Unit : m,
    numBarsSec_Out[i] = Number_Reinforcement_Out[i] = int(
        round(Circumference_Section_Out[i] / Rebar_spacing_VO[i]))  # Number
    numBarsSec_In[i] = Number_Reinforcement_In[i] = int(
        round(Circumference_Section_In[i] / Rebar_spacing_VI[i]))  # Number
    # print('numBarsSec_Out[i]', numBarsSec_Out[i])
    # print('numBarsSec_In[i]', numBarsSec_In[i])
    # print(Rebar_spacing_VO[i])
    # print('Number_Reinforcement ', i, Number_Reinforcement[i], int(round(Number_Reinforcement[i])))
    barAreaSec_Out = Area_Reinforcement_Out = math.pi / 4 * Rebar_diameter_VO ** 2  # m2
    barAreaSec_In = Area_Reinforcement_In = math.pi / 4 * Rebar_diameter_VI ** 2  # m
    # print(barAreaSec_Out[i])
    # print(barAreaSec_In[i])
    # print('Area_Reinforcement ', i, Area_Reinforcement[i])
    ColSecTag[i] = Column_Section_Tag[i] = i
    # print('ColSecTag ', i, ColSecTag[i])
    ro[i] = DSec[i] / 2  # overall (outer) radius of the section
    ri[i] = ro[i] - TSec[i]  # inner radius of the section, only for hollow sections
    # print('ro[', i, ']',DSec[i],'/2 =', ro[i])
    # print('DSec[', i, ']',DSec[i],'/2 - TSec[', i, ']',TSec[i],' =', ri[i])
    # print('ro[i]', ro[i])
    # print('ri[i]', ri[i])
    nfCoreR = 8  # number of radial divisions in the core (number of "rings")
    nfCoverR = 4  # number of radial divisions in the cover
    nfCoreT[i] = int(numBarsSec_Out[i])  # number of theta divisions in the core (number of "wedges")
    # print('nfCoreT[i]', nfCoreT[i])
    nfCoverT[i] = int(numBarsSec_Out[i])  # number of theta divisions in the cover
    # print('nfCoverT[i]', nfCoverT[i])
    rc[i] = ro[i] - CoverSec  # Core radius
    rc2[i] = ri[i] + CoverSec
    # print('rc[i]', rc[i])
    # print('rc2[i]', rc2[i])
    theta_Out[i] = 360.0 / numBarsSec_Out[i]  # Determine angle increment between bars
    theta_In[i] = 360.0 / numBarsSec_In[i]  # Determine angle increment between bars
    # print('theta_Out[i]', theta_Out[i])
    # print('theta_In[i]', theta_In[i])
    # print(type(nfCoreR), type(nfCoreT[i]), type(nfCoverR), type(nfCoverT[i]))
    # print('rc[i]', rc[i])
    rs[1] = math.sqrt(((math.pi * rc[1] ** 2) - 158641.0614) / math.pi)

# -------------------------------------------------------------------------------------------------------
# ----------------------- Define FIBER SECTION properties -----------------------------------------------
# -------------------------------------------------------------------------------------------------------
# Define the fiber section for Circular section
nfCoverT = nfCoverT.tolist()
nfCoreT = nfCoreT.tolist()
numBarsSec_Out = numBarsSec_Out.tolist()
numBarsSec_In = numBarsSec_In.tolist()
for i in range(1, len(element_length)):
    section('Fiber', i)
    if Nodes_coordinate_Y[i] < 100 * Unite:
        IDconcCover = IDconcCoverC40
        IDconcCore = IDconcCoreC40
    else:
        IDconcCover = IDconcCoverC30
        IDconcCore = IDconcCoreC30
    # print('Nodes_coordinate_Y = ', Nodes_coordinate_Y[i], ', IDconcCore =', IDconcCore)
    patch('circ', IDconcCover, nfCoverT[i], 4, *[0, 0], rc[i], ro[i], *[0, 360])
    layer('circ', IDreinf, numBarsSec_Out[i], barAreaSec_Out[i], *[0, 0], rc[i], *[theta_Out[i], 360])
    patch('circ', IDconcCore, nfCoreT[i], 8, *[0, 0], rc2[i], rc[i], *[0, 360])
    layer('circ', IDreinf, numBarsSec_In[i], barAreaSec_In[i], *[0, 0], rc2[i], *[theta_In[i], 360])
    patch('circ', IDconcCover, nfCoverT[i], 4, *[0, 0], ri[i], rc2[i], *[0, 360])
"""
for i in range(8, 10):
    if Nodes_coordinate_Y[i] < 100 * Unite:
        IDconcCover = IDconcCoverC40
        IDconcCore = IDconcCoreC40
    else:
        IDconcCover = IDconcCoverC30
        IDconcCore = IDconcCoreC30
    fib_sec =  [['section', 'Fiber', i, '-GJ', 1.0e6],
                ['patch', 'circ', IDconcCover, nfCoverT[i], 4, *[0, 0], rc[i], ro[i], *[0, 360]],
                ['layer', 'circ', IDreinf, numBarsSec_Out[i], barAreaSec_Out[i], *[0, 0], rc[i], *[0, 360.0-360/numBarsSec_Out[i]]],
                ['patch', 'circ', IDconcCore, nfCoverT[i], 8, *[0, 0], rc2[i], rc[i], *[0, 360]],
                ['layer', 'circ', IDreinf, numBarsSec_In[i], barAreaSec_In[i], *[0, 0], rc2[i], *[0, 360.0-360/numBarsSec_In[i]]],
                ['patch', 'circ', IDconcCover, nfCoverT[i], 4, *[0, 0], ri[i], rc2[i], *[0, 360]]]
    print(fib_sec)
    matcolor = ['r', 'lightgrey', 'gold', 'w', 'w', 'w']
    opsv.plot_fiber_section(fib_sec, matcolor=matcolor)
    plt.axis('equal')
    #plt.savefig('fibsec_rc.png')
    plt.show()
"""
# -------------------------------------------------------------------------------------------------------
# ----------------------- Define NODES properties & BOUNDARY --------------------------------------------
# -------------------------------------------------------------------------------------------------------
for i in range(0, len(element_length)):
    node(i + 1, 0, int(Nodes_coordinate_Y[i]))
    # print('node', i + 1, 0, int(Nodes_coordinate_Y[i]))

# Single point constraints -- Boundary Conditions
fix(1, 1, 1, 1)  # node DX DY RZ

# -------------------------------------------------------------------------------------------------------
# ----------------------- Define ELEMENTS properties ----------------------------------------------------
# -------------------------------------------------------------------------------------------------------
# define geometric transformation: performs a linear geometric transformation of beam
ColTransfTag = 1
geomTransf('PDelta', ColTransfTag)  # associate a tag to transformation

area = np.zeros(len(element_length))
Iz = np.zeros(len(element_length))
Stifness = np.zeros(len(element_length))
for i in range(1, len(element_length)):
    # print('nonlinearBeamColumn', i, *[i, i+1], 5, i, 1)
    # element('nonlinearBeamColumn', eleTag, *eleNodes, numIntgrPts, secTag, transfTag)
    area[i] = math.pi * ((ro[i] ** 2) - (ri[i] ** 2))
    Iz[i] = (math.pi / 4) * ((ro[i] ** 4) - (ri[i] ** 4))
    Stifness[i] = area[i] * EcC30 / element_length[i]
    beamIntegration('Lobatto', i, i, 5)
    ##element('forceBeamColumn', i, *[i, i + 1], ColTransfTag, i)
    ##eelement('nonlinearBeamColumn', i, *[i, i + 1], 5, i, ColTransfTag)
    if Nodes_coordinate_Y[i] < 100 * Unite:
        IDconc = EcC40
    else:
        IDconc = EcC30
    element('elasticBeamColumn', i, *[i, i + 1], area[i], IDconc, Iz[i], ColTransfTag)

# -------------------------------------------------------------------------------------------------------
# ----------------------- Define MASS properties --------------------------------------------------------
# -------------------------------------------------------------------------------------------------------
Density = 2.500 * 1e-9  # ---------- Concrete Density * 1e-9 > t/mm3
Total_Vol = 0
Total_Mass = 0
Mass_Sec = np.zeros(len(element_length))
Vol_Sec = np.zeros(len(element_length))
for i in range(1, len(element_length)):
    Vol_Sec[i] = (math.pi * (ro[i] ** 2 - ri[i] ** 2) * element_length[i])
    Mass_Sec[i] = Vol_Sec[i] * Density
    # print('Volume section', i, Vol_Sec[i])
    # ----------------------------------
    # Additional Mass due to Steel Tubes
    # ----------------------------------
    if i > 14:
        Mass_Sec[i] += 3454 / 10
    else:
        Mass_Sec[i] += 0
    Total_Vol += Vol_Sec[i]
    Total_Mass += Mass_Sec[i]
    mass(i, Mass_Sec[i], Mass_Sec[i], Mass_Sec[i])
    out_mass = nodeMass(i)
    # print(out_mass)
    # print('Mass section', i, Mass_Sec[i])
    # print('mass', i, Mass_Sec[i],  Mass_Sec[i],  Mass_Sec[i])
    # print('Vol_Sec section', i, Vol_Sec[i])
# print('Total Vol of the structure', Total_Vol)
# print('Total Mass of the structure', Total_Mass)

# -------------------------------------------------------------------------------------------------------
# ----------------------- Pulse methode ---------------------------------------------------
# -------------------------------------------------------------------------------------------------------
# """
c = 0.7
kv = 0.65 * 0.05  # 65% alpha(max)
ni = 4 * (1 + c) * kv
Ge = Total_Mass * 9964.0164472819
Gie = np.zeros(len(element_length))
Height = np.zeros(len(element_length))
Fivk = np.zeros(len(element_length))
Sum_Mass_Sec = 0
for i in range(1, len(element_length)):
    Sum_Mass_Sec += Mass_Sec[i]
    Gie[i] = (Total_Mass - Sum_Mass_Sec) * 9964.0164472819
    # print(i)
    # print(Gie[i])
    Fivk[i] = ni * (Gie[i] - (Gie[i] ** 2) / Ge)
    Fivk[0] = 0.75 * Ge * 0.04 * 0.65
    Height[i] = Nodes_coordinate_Y[i] / Unite
    df2 = pd.DataFrame({'Force (N)': Fivk, 'Elevation (m)': Height})
    df2.to_csv('pulse.csv', index=False)
# """

# recorder('Node', '-file', 'GM\\Node_Reaction_Y2_Gravity.out', '-time', '-node', 1, '-dof', 2, 'reaction')
# -------------------------------------------------------------------------------------------------------
# ----------------------- SAVE ODB MODEL FOR PLOTTING ---------------------------------------------------
# -------------------------------------------------------------------------------------------------------
ModelName = 'Chimney'
LoadCaseName = 'Transient'

# opsplt.createODB(ModelName, LoadCaseName, deltaT=0.02, Nmodes=10)
# recorder('Node', '-file', 'Node_Force_Reaction_Y2.out', '-time', '-nodeRange', 1, 25, '-dof', 2, 'reaction')
# recorder('Element', '-file', 'Element_LocalForce_Reaction_Y1.out', '-time', '-eleRange', 1, 24, '-dof', 1, 'localForces')
# -------------------------------------------------------------------------------------------------------
# ----------------------- Define GRAVITY LOAD -----------------------------------------------------------
# -------------------------------------------------------------------------------------------------------
# """
timeSeries('Linear', 400)
pattern('Plain', 400, 400, )
# load(25, 0, -Total_Mass*9806.65, 0)			    # #    nd,  FX,  FY, MZ --  superstructure-weight
for i in range(0, len(element_length) - 1):
    load(i + 1, 0, -Mass_Sec[i + 1] * 9806.65, 0)  # #    nd,  FX,  FY, MZ --  superstructure-weight* 9806.65

constraints('Plain')  # how it handles boundary conditions
numberer('Plain')  # renumber dof's to minimize band-width (optimization), if you want to
system('BandGeneral')  # how to store and solve the system of equations in the analysis
algorithm('Linear')  # use Linear algorithm for linear analysis
integrator('LoadControl', 0.1)  # determine the next time step for an analysis, # apply gravity in 10 steps
analysis('Static')  # define type of analysis static or transient
analyze(10)  # perform gravity analysis
loadConst('-time', 0.0)  # hold gravity constant and restart time

# -------------------------------------------------------------------------------------------------------
# ----------------------- Define damping  ----------------------------------------------------------
# -------------------------------------------------------------------------------------------------------
# Perform an eigenvalue analysis

xDamp = 0.05  # 5% damping ratio
numEigen = 9  # Number of eigenvalue modes
mode_i = 2  # mode 1
mode_j = 5  # mode 3

Lambda = eigen('-fullGenLapack', numEigen)  # eigenvalue modes Lambda
Omega = np.zeros(numEigen)
# -------------------- Elemental Rayleigh Damping
for i in range(0, numEigen):
    Omega[i] = math.pow(Lambda[i], 0.5)
MpropSwitch = 1
KcurrSwitch = 0
KinitSwitch = 0
KcommSwitch = 1
# ----------------------------------------
alphaM = MpropSwitch * xDamp * (2 * Omega[mode_i] * Omega[mode_j]) / (
        Omega[mode_i] + Omega[mode_j])  # M-prop. damping; D = alphaM*M
betaKcurr = KcurrSwitch * 2 * xDamp / (
        Omega[mode_i] + Omega[mode_j])  # K-proportional damping;      +beatKcurr*KCurrent
betaKinit = KinitSwitch * 2 * xDamp / (
        Omega[mode_i] + Omega[mode_j])  # initial-stiffness proportional damping      +beatKinit*Kini
betaKcomm = KcommSwitch * 2 * xDamp / (Omega[mode_i] + Omega[mode_j])  # last-committed K;   +betaKcomm*KlastCommitt

#       (alpha_m, beta_k, beta_k_init, beta_k_comm) # RAYLEIGH  D = αM∗M + βK∗Kcurr + βKinit∗Kinit + βKcomm∗Kcommit
rayleigh(alphaM, betaKcurr, betaKinit, betaKcomm)
# rayleigh(0.0, 0.0, 0.0, 0.000625)

# -------------------------------------------------------------------------------------------------------
# ----------------------- Define GROUND Motion ----------------------------------------------------------
# -------------------------------------------------------------------------------------------------------
# Set some parameters
# ------------------------ SELECTED GM ------------------------->   0,   1,   2
GM_sele = ['waves\GM\RSN1594_CHICHI_TTN051-E', 'waves\GM\RSN6896_DARFIELD_DORCN20W', 'waves\GM\ARTIFICIAL']
scal_fact = [0.4624, 0.1916, 1]  # Scaling Original PEER GM to Intensity 6 design spectrum
"""
# ----------------------------------- PGA Scaling ----------------------------------->  
for i in range(0, len(GM_sele)):
    print('>>>>>> Processing PGA scaling >>>>>> ', GM_sele[i])
    desc, npts, dt, time, inp_acc = processNGAfile(GM_sele[i]+'.AT2', scal_fact[i])
    PGA = np.amax(inp_acc)
    gal_to_g = 0.0010197162129779282
    +----------+--------------+-------------+-------------+
    |          | Intensity 6  | Intensity 7 | Intensity 8 |
    +----------+--------------+-------------+-------------+
    | Minor    | 18(0.04)     | 35(0.08)    | 70(0.16)    |
    | Frequent | 52(0.12)     | 96(0.22)    | 183(0.42)   |
    | Major    | 125(0.28)    | 218(0.50)   | 392(0.90)   |
    +----------+--------------+-------------+-------------+
    PGA_scal = 18 * gal_to_g / PGA          # 18 gal scaling   (6 Minor)
    PGA_scal = 22.5 * gal_to_g / PGA        # 22.5 gal scaling  (6 Design)
    PGA_scal = 125 * gal_to_g / PGA         # 125 gal scaling  (6 Major)
    PGA_scal = 220 * gal_to_g / PGA         # 220 gal scaling  (7 Major)
    print('Max element from Numpy Array : ', PGA_scal)
"""
scal_fact_PGA = [1.37157503147026, 1.3747569775679702, 1.2916843771632733]  # Scaling PGA Intensity 6 to Various
# Active_GM = np.zeros(len(GM_sele))


iter = 0
# i = 2
gal = [18, 52, 125, 218, 392]
for i in range(0, len(GM_sele)):
    for j in range(0, len(gal)):
        Active_GM = i
        record = GM_sele[i]
        # Perform the conversion from SMD record to OpenSees record
        dt, nPts = ReadRecord.ReadRecord(record + '.at2', record + '.dat')
        print((record + '.at2', record + '.dat'))
        # Uniform EXCITATION: acceleration input
        GMfact = scal_fact[i] * scal_fact_PGA[i] * 9806.65 * (
                    gal[j] / 18)  # data in input file is in g Units -- Scale factor + (9806.65) ---> g to mm/s2
        # GMdirection = [1, 2]  #  1=UX, 2=UY, 3=UZ, 4=RX, 5=RY, 6=RZ
        # GM_direction_fact = [1, 0.65]      # Horizontal 1 = 100%, Vertical 0.65 = 65%

        # -------------------------------------------------------------------------------------------------------
        # ----------------------- Files Recorder  ---------------------------------------------------------------
        # -------------------------------------------------------------------------------------------------------
        if gal[j] == 18:
            nfolder = '6_Minor'
        elif gal[j] == 52:
            nfolder = '6_Frequent'
        elif gal[j] == 125:
            nfolder = '6_Major'
        elif gal[j] == 218:
            nfolder = '7_Major'
        elif gal[j] == 392:
            nfolder = '8_Major'
        else:
            pass

        if Active_GM == 0:  ############### GM 1 Recorder ######################
            print('Active_GM1')
            recorder('Node', '-file', f'GM\\{nfolder}\\Node_Disp_Y2_GM1.out', '-time', '-nodeRange', 1, 25, '-dof', 2,
                     'disp')
            recorder('Node', '-file', f'GM\\{nfolder}\\Node_Acc_Reaction_Y2_GM1.out', '-time', '-nodeRange', 1, 25,
                     '-dof', 2, 'accel')
        elif Active_GM == 1:  ############### GM 2 Recorder ######################
            print('Active_GM2')
            recorder('Node', '-file', f'GM\\{nfolder}\\Node_Disp_Y2_GM2.out', '-time', '-nodeRange', 1, 25, '-dof', 2,
                     'disp')
            recorder('Node', '-file', f'GM\\{nfolder}\\Node_Acc_Reaction_Y2_GM2.out', '-time', '-nodeRange', 1, 25,
                     '-dof', 2, 'accel')
        elif Active_GM == 2:  ############### GM 3 Recorder ######################
            print('Active_GM3')
            recorder('Node', '-file', f'GM\\{nfolder}\\Node_Disp_Y2_GM3.out', '-time', '-nodeRange', 1, 25, '-dof', 2,
                     'disp')
            recorder('Node', '-file', f'GM\\{nfolder}\\Node_Acc_Reaction_Y2_GM3.out', '-time', '-nodeRange', 1, 25,
                     '-dof', 2, 'accel')
        else:
            pass
        # -------------------------------------------------------------------------------------------------------
        # ----------------------- DYNAMIC ground-motion analysis (Transient analysis) ---------------------------
        # -------------------------------------------------------------------------------------------------------
        # ------------------ Horizontal Direction ------------------------------------------------------
        # timeSeries('Path', i+1, '-dt', dt, '-filePath', record + '.dat', '-factor', GMfact * 1)
        # pattern('UniformExcitation', 2, 1, '-accel', 2)  # 1=UX, 2=UY, 3=UZ, 4=RX, 5=RY, 6=RZ
        # ------------------ Vertical Direction ------------------------------------------------------
        print('iter--------------------', iter)
        timeSeries('Path', iter, '-dt', dt, '-filePath', record + '.dat', '-factor', GMfact * 0.65)
        print('Path', iter, '-dt', dt, '-filePath', record + '.dat', '-factor', GMfact * 0.65)
        pattern('UniformExcitation', iter, 2, '-accel',
                iter)  # 1=UX, 2=UY, 3=UZ, 4=RX, 5=RY, 6=RZ pattern('UniformExcitation', patternTag, dir, '-disp', dispSeriesTag, '-vel', velSeriesTag, '-accel', accelSeriesTag, '-vel0', vel0, '-fact', fact)
        print('UniformExcitation', iter, 2, '-accel', iter)
        # Create the transient analysis------------------------------------
        wipeAnalysis()
        constraints('Transformation')
        numberer('Plain')
        system('BandGeneral')
        test('EnergyIncr', 1e-6, 200, 1)
        algorithm('ModifiedNewton')
        integrator('Newmark', 0.5, 0.25)
        analysis('Transient')
        print(nPts)
        print(dt)
        iter += 1
        analyze(nPts, dt)




# -------------------------------------------------------------------------------------------------------
# ----------------------- Response Spectrum Method ----------------------------------------------------
# -------------------------------------------------------------------------------------------------------

# from Modal import ModalAnalysis
# ModalAnalysis(5, pflag=1, outname='modalAnalysis')
# printA('-file', 'matrix_transient.out')

"""
ok = analyze(nPts, dt)

Tol = 1e-8  # convergence tolerance for test
te = {1: 'NormDispIncr', 2: 'RelativeEnergyIncr', 4: 'RelativeNormUnbalance', 5: 'RelativeNormDispIncr',
      6: 'NormUnbalance'}
algo = {1: 'KrylovNewton', 2: 'SecantNewton', 4: 'RaphsonNewton', 5: 'PeriodicNewton', 6: 'BFGS', 7: 'Broyden',
        8: 'NewtonLineSearch'}
for i in te:
    for j in algo:
        if ok != 0:
            if j < 4:
                algorithm(algo[j], '-initial')

            else:
                algorithm(algo[j], '-initial')

            test(te[i], Tol, 200)
            ok = analyze(nPts, dt)
            print(te[i], algo[j], ok)
            if ok == 0:
                break
        else:
            continue
#"""

# -------------------------------------------------------------------------------------------------------
# ----------------------- Plot Model and Deformation ----------------------------------------------------
# -------------------------------------------------------------------------------------------------------

# plt.plot(reaction, Time)
# plt.title(filename)
# plt.grid(linestyle='--', linewidth=0.5)
# plt.ylabel('Displacement H (cm/s2)')
# plt.xlabel('Time (s)')
# plt.show()

# plt.figure()
# opsplt.plot_model("nodes", "elements")
# opsplt.plot_modeshape(1, 100)
# opsplt.plot_modeshape(2, 100)
# opsplt.plot_modeshape(3, 100)
# opsplt.plot_deformedshape(Model=ModelName, LoadCase=LoadCaseName, scale=100, overlap="yes")  # tstep=100.0,
# plt.show()
# opsv.plot_mode_shape(3, sfac=200000000, unDefoFlag=0, endDispFlag=0)
# ani = opsplt.animate_deformedshape(ModelName, LoadCaseName, dt=dt)
# plt.show()

# opsplt.plot_modeshape(2, 200, Model = ModelName)
# opsplt.plot_deformedshape(Model = ModelName, LoadCase = LoadCaseName)
# anime = opsplt.animate_deformedshape(Model = ModelName, LoadCase = LoadCaseName, dt = dt, tStart=10.0, tEnd=20.0, scale=200, Movie="Dynamic")



