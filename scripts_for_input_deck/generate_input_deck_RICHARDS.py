import sys

nx = 1; dx = 1; lx = 1
ny = 1; dy = 1; ly = 1
nz = 1; dz = 1; lz = 1
filename = ''
for k in range(len(sys.argv)-1):
  k = k + 1
  if sys.argv[k] == '-nx':
    nx = int(float(sys.argv[k+1]))
  if sys.argv[k] == '-ny':
    ny = int(float(sys.argv[k+1]))
  if sys.argv[k] == '-nz':
    nz = int(float(sys.argv[k+1]))
  if sys.argv[k] == '-dx':
    dx = float(sys.argv[k+1])
  if sys.argv[k] == '-dy':
    dy = float(sys.argv[k+1])
  if sys.argv[k] == '-dz':
    dz = float(sys.argv[k+1])
  if sys.argv[k] == '-lx':
    lx = float(sys.argv[k+1])
  if sys.argv[k] == '-ly':
    ly = float(sys.argv[k+1])
  if sys.argv[k] == '-lz':
    lz = float(sys.argv[k+1])
  if sys.argv[k] == '-input_prefix':
    filename = sys.argv[k+1] + '.in'
    
f = open(filename,'w')

simulation = """ 
SIMULATION
  SIMULATION_TYPE SUBSURFACE
  PROCESS_MODELS
    SUBSURFACE_FLOW flow
      MODE RICHARDS
    /
  /
END

"""
grid = """
GRID
  TYPE structured
  NXYZ """ + str(nx) + " " + str(ny) + " " + str(nz) + """
  DXYZ
   """ + str(dx) + """d0
   """ + str(dy) + """d0
   """ + str(dz) + """d0
  END
END
"""
region = """
REGION all
  COORDINATES
    0.0d0 0.0d0 0.0d0
    """ + str(lx) + "d0 " + str(ly) + "d0 " + str(lz) + """d0
  /
END

REGION left
  COORDINATES
    0.d0   0.d0   0.d0
    """ + str(lx*0.40) + "d0 " + str(ly) + "d0 " + str(lz) + """d0
  /
END

REGION right
  COORDINATES
    """ + str(lx*0.40) + """d0 0.0d0 0.0d0
    """ + str(lx) + "d0 " + str(ly) + "d0 " + str(lz) + """d0
  /
END

REGION left_face
  FACE WEST
  COORDINATES
    0.d0 0.d0 0.d0
    0.d0 """ + str(ly) + "d0 " + str(lz) + """d0
  /
END

REGION right_face
  FACE EAST
  COORDINATES
    """ + str(lx) + """d0 0.d0 0.d0
    """ + str(lx) + "d0 " + str(ly) + "d0 " + str(lz) + """d0
  /
END
"""
mat_prop = """
MATERIAL_PROPERTY beam_left
  ID 1
  POROSITY 0.5
  TORTUOSITY 1.d0
  ROCK_DENSITY 2.8E3 kg/m^3
  CHARACTERISTIC_CURVES cc1
  PERMEABILITY
    PERM_X 1.d-12
    PERM_Y 1.d-12
    PERM_Z 1.d-12
  /
END

MATERIAL_PROPERTY beam_right
  ID 2
  POROSITY 0.5
  TORTUOSITY 1.d0
  ROCK_DENSITY 2.8E3 kg/m^3
  CHARACTERISTIC_CURVES cc1
  PERMEABILITY
    PERM_X 3.d-12
    PERM_Y 3.d-12
    PERM_Z 3.d-12
  /
END
"""
char_curves = """
CHARACTERISTIC_CURVES cc1
  SATURATION_FUNCTION VAN_GENUCHTEN
    LIQUID_RESIDUAL_SATURATION 0.5d-1
    M 0.75
    ALPHA 1.d-3
  /
  PERMEABILITY_FUNCTION MUALEM
    LIQUID_RESIDUAL_SATURATION 0.1d0
    M 0.5d0
  /
  PERMEABILITY_FUNCTION MUALEM_VG_GAS
    LIQUID_RESIDUAL_SATURATION 0.1d0
    GAS_RESIDUAL_SATURATION 0.1d0
    M 0.5d0
  /
END
"""
strata = """
STRATA
  REGION left
  MATERIAL beam_left
END

STRATA
  REGION right
  MATERIAL beam_right
END
"""
fluid_prop = """
FLUID_PROPERTY
  DIFFUSION_COEFFICIENT 1.d-9
END
"""
eoswater = """
EOS WATER
 DENSITY CONSTANT 1000.d0 kg/m^3
 VISCOSITY CONSTANT 1.0d-3 Pa-s
END
"""
time = """
TIME
  FINAL_TIME 10 day
  INITIAL_TIMESTEP_SIZE 1.d-4 day
  MAXIMUM_TIMESTEP_SIZE 0.25 day at 0.d0 day
END
"""
output = """
OUTPUT
  SNAPSHOT_FILE
    NO_PRINT_INITIAL
    #NO_PRINT_FINAL
    FORMAT HDF5
  /
END
"""
flow_cond = """
FLOW_CONDITION initial
  TYPE
    PRESSURE dirichlet
  /
  PRESSURE 1.0 MPa
END

FLOW_CONDITION left_face
  TYPE
    PRESSURE dirichlet
  /
  PRESSURE 1.0 MPa
END

FLOW_CONDITION right_face
  TYPE
    FLUX neumann
  /
  FLUX +1.50d-5 m/sec
END
"""
init_cond = """
INITIAL_CONDITION initial_left
  REGION left
  FLOW_CONDITION initial
END

INITIAL_CONDITION initial_right
  REGION right
  FLOW_CONDITION initial
END

BOUNDARY_CONDITION left_face
  REGION left_face
  FLOW_CONDITION left_face
END

BOUNDARY_CONDITION right_face
  REGION right_face
  FLOW_CONDITION right_face
END
"""
  
  
f.write(simulation)
f.write("SUBSURFACE\n")
f.write(grid)
f.write(region)
f.write(mat_prop)
f.write(char_curves)
f.write(strata)
f.write(fluid_prop)
f.write(eoswater)
f.write(time)
f.write(output)
f.write(flow_cond)
f.write(init_cond)
f.write("\nEND_SUBSURFACE\n") 
# Must have \n at the end, or HDF5 output will not work
