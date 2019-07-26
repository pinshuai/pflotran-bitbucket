# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 14:44:16 2019

@author: gehammo
"""

import sys
import shutil
import os
import re

filename = 'pflotran.in'
f = open(filename,'w')

nx = 200
ny = 1
nz = 1
length_x = 100.
length_y = 1.
length_z = 1.
velx = 1.
vely = 0.
velz = 5.

simulation_block = '''
SIMULATION
  SIMULATION_TYPE SUBSURFACE
  PROCESS_MODELS
    SUBSURFACE_TRANSPORT transport
      GLOBAL_IMPLICIT
    /
  /
END

SUBSURFACE
'''

velocity_block = '''
SPECIFIED_VELOCITY
  UNIFORM? YES
  DATASET {} {} {} m/yr
END
'''


grid_template = '''
GRID
  TYPE STRUCTURED
  NXYZ {} {} {}
  BOUNDS
    0.d0 0.d0 0.d0
    {} {} {}
  /
END
'''

f.write(simulation_block)
f.write(velocity_block.format(velx,vely,velz))
f.write(grid_template.format(nx,ny,nz,length_x,length_y,length_z))
f.write('END_SUBSURFACE')

print('done')
f.close()