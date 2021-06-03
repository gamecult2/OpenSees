print("============================================")
print("Start Chimney EQ ground motion with gravity ")

import openseespy.postprocessing.ops_vis as opsv
from openseespy.opensees import *
import numpy as np
import matplotlib.pyplot as plt
import math

# SET UP -----------------------------------------------------
wipe()  # clear opensees model
model('basic', '-ndm', 2, '-ndf', 3)  # 2 dimensions, 3 dof per node

# ------------------------------------------------------------
# ----------------------- CONCRETE C40 -----------------------
# ------------------------------------------------------------
IDconcCoreC40 = 1                                    # material ID tag -- confined core concrete
IDconcCoverC40 = 2                                   # material ID tag -- unconfined cover concrete
# nominal concrete compressive strength
fcC40 = -34.99                                       # CONCRETE Compressive Strength, MPa   (+Tension, -Compression)
EcC40 = 31500                                        # Concrete Elastic Modulus
# confined concrete
KfcC40 = 1.3                                         # ratio of confined to unconfined concrete strength
fc1CC40 = KfcC40 * fcC40                             # CONFINED concrete (mander model), maximum stress
eps1CC40 = 2 * fc1CC40 / EcC40                       # strain at maximum stress
fc2CC40 = 0.2 * fc1CC40                              # ultimate stress
eps2CC40 = 5 * eps1CC40                              # strain at ultimate stress
# unconfined concrete
fc1UC40 = fcC40                                      # UNCONFINED concrete (todeschini parabolic model), maximum stress
eps1UC40 = -0.003                                    # strain at maximum strength of unconfined concrete
fc2UC40 = 0.2 * fc1UC40                              # ultimate stress
eps2UC40 = -0.01                                     # strain at ultimate stress
lamC40 = 0.1                                         # ratio between unloading slope at $eps2 and initial slope $Ec
# tensile-strength properties
ftCC40 = -0.14 * fc1CC40                             # tensile strength +tension
ftUC40 = -0.14 * fc1UC40                             # tensile strength +tension
EtsC40 = ftUC40 / 0.002                              # tension softening stiffness

# ------------------------------------------------------------
# ----------------------- CONCRETE C30 -----------------------
# ------------------------------------------------------------
IDconcCoreC30 = 3  # material ID tag -- confined core concrete
IDconcCoverC30 = 4  # material ID tag -- unconfined cover concrete
# nominal concrete compressive strength
fcC30 = -34.99  # CONCRETE Compressive Strength, MPa   (+Tension, -Compression)
EcC30 = 31500  # Concrete Elastic Modulus
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
ftUC30 = -0.14 * fc1UC30 # tensile strength +tension
EtsC30 = ftUC30 / 0.002  # tension softening stiffness

# ------------------------------------------------------------
# ----------------------- STEEL ------------------------------
# ------------------------------------------------------------
IDreinf = 5  # material ID tag -- reinforcement

Fy = 413.7  # STEEL yield stress
Es = 199947.9  # modulus of steel
Bs = 0.01  # strain-hardening ratio
R0 = 18.5  # control the transition from elastic to plastic branches
cR1 = 0.925  # control the transition from elastic to plastic branches
cR2 = 0.15	# control the transition from elastic to plastic branches

# Define materials for nonlinear columns
# -------------------------------------------------------------------
# CORE CONCRETE C40 (confined)
uniaxialMaterial('Concrete02', IDconcCoreC40, fc1CC40, eps1CC40, fc2CC40, eps2CC40, lamC40, ftCC40, EtsC40)
# COVER CONCRETE C40 (unconfined)
uniaxialMaterial('Concrete02', IDconcCoverC40, fc1CC40, eps1CC40, fc2CC40, eps2CC40, lamC40, ftCC40, EtsC40)
# ------------------------------------------------------------------
# CORE CONCRETE C30 (confined)
uniaxialMaterial('Concrete02', IDconcCoreC30, fc1CC30, eps1CC30, fc2CC30, eps2CC30, lamC30, ftCC30, EtsC30)
# COVER CONCRETE C30 (unconfined)
uniaxialMaterial('Concrete02', IDconcCoverC30, fc1CC30, eps1CC30, fc2CC30, eps2CC30, lamC30, ftCC30, EtsC30)
# ------------------------------------------------------------------
# STEEL Reinforcing steel
uniaxialMaterial('Steel02', IDreinf, Fy, Es, Bs, R0, cR1, cR2)
# END MATERIAL parameters ------------------------------------------


# ------------------------------------------------------------------
# ----------------------- Different GEOMETRY Coordinate-------------
# ------------------------------------------------------------------
#                          1   2   3   4   5   6   7   8   9   10  11  12  13  14  15  16  17  18  19  20  21  22  23  24
element_length = np.array([10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 4])  # Number of element 25
#                              1  2   3   4   5   6   7   8   9   10  11   12   13   14   15   16   17   18   19   20   21   22   23   24   25
Nodes_coordinate_Y = np.array([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 234])
# Section Diameter ---------------         1     2     3     4     5     6     7     8     9     10    11    12    13    14    15    16    17    18    19    20    21    22    23    24
DSec = Column_Section_Diameter = np.array([27.3, 26.7, 26.1, 25.5, 24.9, 24.4, 23.7, 23.5, 23.1, 22.7, 22.3, 22.1, 21.9, 21.7, 21.5, 21.3, 21.1, 20.9, 20.9, 20.9, 20.9, 20.9, 20.9, 20.9])
# Section Thickness ---------------         1     2     3     4     5     6     7     8     9     10    11    12    13    14    15    16    17    18    19    20    21    22    23    24
TSec = Column_Section_Thickness = np.array([0.75, 0.73, 0.71, 0.69, 0.67, 0.65, 0.65, 0.65, 0.63, 0.60, 0.57, 0.55, 0.53, 0.51, 0.49, 0.45, 0.43, 0.41, 0.39, 0.39, 0.39, 0.35, 0.35, 0.35])

# Vertical reinforcement on the outer side of the tube wall (HRB400 grade is adopted for the reinforcement)
Rebar_spacing_VO = np.array([0.205, 0.205, 0.205, 0.215, 0.215, 0.225, 0.225, 0.225, 0.225, 0.235, 0.235, 0.235, 0.235, 0.235, 0.235, 0.235, 0.235, 0.235, 0.235, 0.235, 0.235, 0.235, 0.235, 0.235, 0.235, 0.235, 0.235, 0.235])
Rebar_diameter_VO = np.array([0.022, 0.022, 0.022, 0.022, 0.022, 0.022, 0.022, 0.022, 0.022, 0.022, 0.022, 0.022, 0.022, 0.022, 0.022, 0.022, 0.020, 0.020, 0.020, 0.020, 0.020, 0.018, 0.018, 0.018, 0.018, 0.018, 0.018, 0.018])

# Vertical reinforcement on the inner side of the tube wall (HRB400 grade is adopted for the reinforcement)
Rebar_spacing_VI = np.full(len(element_length), 0.200)
Rebar_diameter_VI = np.array(
    [0.022, 0.022, 0.022, 0.022, 0.022, 0.022, 0.022, 0.022, 0.022, 0.022, 0.022, 0.022, 0.022, 0.022, 0.022, 0.022,
     0.020, 0.020, 0.020, 0.020, 0.020, 0.018, 0.018, 0.018, 0.018, 0.018, 0.018, 0.018])

# Circumferential reinforcement on the outside of the wall (HRB400 grade is adopted for the reinforcement)
#Rebar_spacing_CO = np.full(len(element_length), 0.200)
#Rebar_diameter_CO = np.array(
#    [0.022, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.018, 0.018,
#    0.018, 0.018, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016])

# Circumferential steel bars on the inner side of the tube wall (HRB400 grade is adopted for the reinforcement)
#Rebar_spacing_CI = np.array(
#    [0.200, 0.200, 0.200, 0.200, 0.200, 0.200, 0.200, 0.200, 0.200, 0.200, 0.200, 0.200, 0.200, 0.200, 0.200, 0.200,
#     0.200, 0.200, 0.200, 0.200, 0.200, 0.200, 0.190, 0.175, 0.165, 0.140, 0.135, 0.135])
#Rebar_diameter_CI = np.array(
#    [0.018, 0.018, 0.018, 0.018, 0.018, 0.018, 0.018, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.016,
#     0.016, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014, 0.014])

Circumference_Section = np.zeros(len(element_length))
numBarsSec = Number_Reinforcement = np.zeros(len(element_length), dtype='i')
barAreaSec = Area_Reinforcement = np.zeros(len(element_length))
ColSecTag = Column_Section_Tag = np.zeros(len(element_length), dtype='i')
ro = np.zeros(len(element_length))
ri = np.zeros(len(element_length))
rc = np.zeros(len(element_length))
rc2 = np.zeros(len(element_length))
theta = np.zeros(len(element_length))
nfCoreT = np.zeros(len(element_length), dtype='i')
nfCoverT = np.zeros(len(element_length), dtype='i')

CoverSec = 0.05
for i in range(1, len(element_length)):
    Circumference_Section[i] = DSec[i] * math.pi  # m
    #print('Element ', i, Circumference_Section[i])
    numBarsSec[i] = Number_Reinforcement[i] = int(round(Circumference_Section[i] / Rebar_spacing_VO[i]))  # Number
    #print(Rebar_spacing_VO[i])
    #print('Number_Reinforcement ', i, Number_Reinforcement[i], int(round(Number_Reinforcement[i])))
    barAreaSec = Area_Reinforcement = math.pi / 4 * Rebar_diameter_VO ** 2  # m2
    # print('Area_Reinforcement ', i, Area_Reinforcement[i])
    ColSecTag[i] = Column_Section_Tag[i] = i
    # print('ColSecTag ', i, ColSecTag[i])
    ro[i] = DSec[i] / 2  # overall (outer) radius of the section
    ri[i] = ro[i] - TSec[i]  # inner radius of the section, only for hollow sections
    #print('ro[', i, ']',DSec[i],'/2 =', ro[i])
    #print('DSec[', i, ']',DSec[i],'/2 - TSec[', i, ']',TSec[i],' =', ri[i])
    #print('ro[i]', ro[i])
    #print('ri[i]', ri[i])
    nfCoreR = 8                             # number of radial divisions in the core (number of "rings")
    nfCoreT[i] = numBarsSec[i]              # number of theta divisions in the core (number of "wedges")
    nfCoverR = 4                            # number of radial divisions in the cover
    nfCoverT[i] = numBarsSec[i]             # number of theta divisions in the cover
    rc[i] = ro[i] - CoverSec                # Core radius
    rc2[i] = ri[i] + CoverSec
    theta[i] = 360.0 / numBarsSec[i]        # Determine angle increment between bars
    #print(type(nfCoreR), type(nfCoreT[i]), type(nfCoverR), type(nfCoverT[i]))

# -----------------------------------------------------------------------------------
# ----------------------- Define FIBER SECTION properties ---------------------------
# -----------------------------------------------------------------------------------
# Define the fiber section for Circular section
for i in range(1, len(element_length)):
    section('Fiber', ColSecTag[i])
    patch('circ', IDconcCoverC40, nfCoverT[i], nfCoverR, 0, 0, rc[i], ro[i], 0, 360)  # Define the cover patch
    layer('circ', IDreinf, numBarsSec[i], barAreaSec[i], 0., 0., rc[i], theta[i], 360)  # Define the reinforcing layer
    patch('circ', IDconcCoreC40, nfCoreT[i], nfCoreR, 0., 0., rc2[i], rc[i], 0., 360)  # Define the core patch
    layer('circ', IDreinf, numBarsSec[i], barAreaSec[i], 0., 0., rc2[i], theta[i], 360)  # Define the reinforcing layer
    patch('circ', IDconcCoverC40, nfCoverT[i], nfCoverR, 0., 0., ri[i], rc2[i], 0., 360)  # Define the cover patch
# END FIBER SECTION properties -------------------------------------------------------

# -----------------------------------------------------------------------------------
# ----------------------- Define GEOMETRY properties --------------------------------
# -----------------------------------------------------------------------------------
for i in range(0, len(element_length)):
    node(i + 1, 0, int(Nodes_coordinate_Y[i]))
    #print('node', i+1, 0, int(Nodes_coordinate_Y[i]))

# Single point constraints -- Boundary Conditions
fix(1, 1, 1, 1)  # node DX DY RZ

# -----------------------------------------------------------------------------------
# ----------------------- Define ELEMENTS properties --------------------------------
# -----------------------------------------------------------------------------------
# define geometric transformation: performs a linear geometric transformation of beam stiffness and resisting force from the basic system to the global-coordinate system
ColTransfTag = 1
geomTransf('PDelta',  ColTransfTag)  		             # associate a tag to transformation

for i in range(0, len(element_length)):
    print('nonlinearBeamColumn', i, *[i, i+1], 5, ColSecTag[i], 1)
    element('nonlinearBeamColumn', i, *[i, i+1], 5, ColSecTag[i], 1)

# fib_sec = [['section', 'Fiber', 1, '-GJ', 1.0e6],
#             ['patch', 'circ', IDconcCoverC40, nfCoverT[i], nfCoverR, 0., 0., rc[i], ro[i], 0., 360],
#             ['layer', 'circ', IDreinf, numBarsSec[i], barAreaSec[i], 0., 0., rc[i], theta[i], 360],
#             ['patch', 'circ', IDconcCoreC40, nfCoreT[i], nfCoreR, 0., 0., rc2[i], rc[i], 0., 360],
#             ['layer', 'circ', IDreinf, numBarsSec[i], barAreaSec[i], 0., 0., rc2[i], theta[i], 360],
#             ['patch', 'circ', IDconcCoverC40, nfCoverT[i], nfCoverR, 0., 0., ri[i], rc2[i], 0., 360]]
# matcolor = ['r', 'lightgrey', 'gold', 'w', 'w', 'w']
# opsv.plot_fiber_section(fib_sec, matcolor=matcolor)
# plt.axis('equal')
# plt.savefig('fibsec_rc.png')
# plt.show()
