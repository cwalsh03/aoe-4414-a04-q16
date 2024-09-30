# ecef_to_sez.py
#
# Usage: python3 ecef_to_sez.py o_x_km o_y_km o_z_km x_km y_km z_km
#  o_ represents the ECEF origin of the SEZ frame, and the other coordinates represent the ECEF position.
# 
# Parameters:
#  arg1: description of argument 1
#  arg2: description of argument 2
#  ...
# Output:
#  A description of the script output
#
# Written by First Last
# Other contributors: None
#
# Optional license statement, e.g., See the LICENSE file for the license.

# import Python modules
# e.g., import math # math module
import sys # argv
import math
import numpy as np

# "constants"
R_E_KM = 6378.137
E_E = 0.081819221456

# helper functions
## calculated denominator
def calc_denom(ecc, lat_rad):
    return math.sqrt(1.0-(E_E**2)*(math.sin(lat_rad)**2)
)

# initialize script arguments
o_x_km = float('nan') # x coord of ECEF origin
o_y_km = float('nan') # y coord of ECEF origin
o_z_km = float('nan') # z coord of ECEF origin
x_km = float('nan') # x coord of ECEF position 
y_km = float('nan') # x coord of ECEF position
z_km= float('nan') # x coord of ECEF position

# parse script arguments
if len(sys.argv)==7:
    o_x_km = float(sys.argv[1]) 
    o_y_km = float(sys.argv[2]) 
    o_z_km = float(sys.argv[3])
    x_km = float(sys.argv[4]) 
    y_km = float(sys.argv[5]) 
    z_km= float(sys.argv[6]) 
else:
  print(\
   'Usage: '\
   'python3 o_x_km o_y_km o_z_km x_km y_km z_km'\
  )
  exit()

# write script below this line

r_o_ECEF = [ [o_x_km],
             [o_y_km],
             [o_z_km]
]

r_ECEF = [ [x_km],
           [y_km],
           [z_km]
]

rECEF = np.subtract(r_ECEF , r_o_ECEF)

r_x_km = rECEF[0][0]
r_y_km = rECEF[1][0] 
r_z_km = rECEF[2][0]

# calculate longitude
lon_rad = math.atan2(o_y_km,o_x_km)
lon_deg = lon_rad*180.0/math.pi

# initialize lat_rad, r_lon_km, r_z_km
lat_rad = math.asin(o_z_km/math.sqrt(o_x_km**2+o_y_km**2+o_z_km**2))
o_lon_km = math.sqrt(o_x_km**2+o_y_km**2)
prev_lat_rad = float('nan')

# iteratively find latitude
c_E = float('nan')
count = 0
while (math.isnan(prev_lat_rad) or abs(lat_rad-prev_lat_rad)>10e-7) and count<5:
  denom = calc_denom(E_E,lat_rad)
  c_E = R_E_KM/denom
  prev_lat_rad = lat_rad
  lat_rad = math.atan((o_z_km+c_E*(E_E**2)*math.sin(lat_rad))/o_lon_km)
  count = count+1
  
denom = calc_denom(E_E, lon_rad)
c_E = R_E_KM / denom
S_E = (R_E_KM*(1 - E_E*E_E)) / denom

# calculate hae
hae_km = o_lon_km/math.cos(lat_rad)-c_E

phi = lon_rad # latitude in rad
theta = lat_rad # longitude in rad

# denom = calc_denom(E_E, phi)
# C_E = R_E_KM / denom
# S_E = (R_E_KM*(1 - E_E*E_E)) / denom

# Rotation matrices - apply latitude rotation first, then longitude
Rzinv = np.array([[math.sin(lat_rad),  0, -math.cos(lat_rad)],
                  [0,                  1,  0],
                  [math.cos(lat_rad),  0,  math.sin(lat_rad)]])

Ryinv = np.array([[math.cos(lon_rad), math.sin(lon_rad), 0],
                  [-math.sin(lon_rad), math.cos(lon_rad), 0],
                  [0,                  0,                1]])

# Multiply matrices
r_SEZ = Rzinv @ Ryinv @ rECEF


s_km = r_SEZ[0][0]
e_km = r_SEZ[1][0]
zz_km = r_SEZ[2][0]


print(s_km)
print(e_km)
print(zz_km)