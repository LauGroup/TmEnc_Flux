#draws a circle around the pore in pymol. run after ionflow.py to verify coordinates.

import math
import pymol
from pymol.cgo import *
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis import distances
from MDAnalysis.lib.log import ProgressBar


# backup constants 
cyl_rad = 7.5 #angstrom
cyl_height = 5.0 #angstrom


print("ENCAPSULIN_FLUX_COUNTER\n\n")
mda.start_logging()
u = mda.Universe("M26_ions_1264.prmtop", "M26_ions_1264_combined_imaged_full.nc")

#select all the F180s
center_pore_ag = u.select_atoms("(resid 180 444 708 972 1236) and name CA")
assert(len(set(center_pore_ag.resnames)) == 1)
assert(list(set(center_pore_ag.resnames))[0] == "PHE")

#compute minimum distance between F180 Calphas  and set the cylinder radius to this value
self_distances = distances.self_distance_array(center_pore_ag.positions)
max_dist = np.max(self_distances)
print(f"Max distance between F180s is {max_dist}")

center_pore_com = center_pore_ag.center_of_mass()
print(f"Pore COM is  {center_pore_com}\nSetting cyl_center to this value\n")
    

def cgoCircle(x, y, z, r=8.0, cr=1.0, cg=0.4, cb=0.8, w=2.0):
  """
  Create a CGO circle
  PARAMS
        x, y, z
          X, Y and Z coordinates of the origin
        r
          Radius of the circle
        cr, cg, cb
          Color triplet, [r,g,b] where r,g,b are all [0.0,1.0].
        w
          Line width of the circle
  RETURNS
        the CGO object (it also loads it into PyMOL, too).
  """
  x = float(x)
  y = float(y)
  z = float(z)
  r = abs(float(r))
  cr = abs(float(cr))
  cg = abs(float(cg))
  cb = abs(float(cb))
  w = float(w)

  obj = [ BEGIN, LINES, COLOR, cr, cg, cb ]
  for i in range(180):
        obj.append( VERTEX )
        obj.append(r*math.cos(i) + x )
        obj.append(r*math.sin(i) + y )
        obj.append(z)
        obj.append( VERTEX )
        obj.append(r*math.cos(i+0.1) + x )
        obj.append(r*math.sin(i+0.1) + y )
        obj.append(z)
  obj.append(END)
 
  cName = cmd.get_unused_name("circle_")
  cmd.load_cgo( obj, cName )
  cmd.set("cgo_line_width", w, cName )
  return obj


r1,g1,b1 = 0,1,0
r2,g2,b2 = 0,1,1
cgoCircle(center_pore_com[0],center_pore_com[1],center_pore_com[2],2*cyl_rad)
cmd.load_cgo( [ ALPHA, 0.7, CYLINDER, center_pore_com[0],center_pore_com[1],center_pore_com[2]-cyl_height, center_pore_com[0],center_pore_com[1],center_pore_com[2]+cyl_height, cyl_rad*2, r1, g1, b1, r2, g2, b2 ], "cylinder" )