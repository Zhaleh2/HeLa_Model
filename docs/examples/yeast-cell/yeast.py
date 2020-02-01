# Import general objects
import lm
import pyLM
import pyLM.units
from pyLM.units import *
import pyLM.RDME
from pyLM.RDME import RDMESimulation
import lmarray

# Import functions just for this simulation
from yeast_geometry import geometryModel
from reactions import reactionModel
from reactions import particleModel
from diffusion import diffusionModel


# Filename
filename =  "Yeast-full"

#######################
# Template Simulation #
#######################
print("Creating Species...")
def initRDME():
    latticeSpacing = 28.7 # nm
    Dfastest = 6.1e-13 # The fastest diffusion coefficint is required for the calculation of the timestep
    Tsim = 15 # Simulation time in seconds
    sim = RDMESimulation(dimensions=micron(5.51, 5.51, 5.51), # Dimensions correspond to a 5.25um x 3.5um x 3.5um ellipsoidal cell
                         spacing=nm(latticeSpacing), 
                         defaultRegion='extra')
    
    # Define simulation timing parameters
    dt = (nm(latticeSpacing)**2)/(2*Dfastest)  
    print('dt =', dt)
    sim.setTimestep(dt)
    sim.setWriteInterval(10*dt) # Write interval of species
    sim.setLatticeWriteInterval(100*dt) # Write interval of the full lattice 
    sim.setSimulationTime(Tsim) 
    
    ###################
    # Specify Species #
    ###################
    print("Defining Species...")
    species = ['gene','premRNA','S','SpremRNA','mRNA']
    sim.defineSpecies(species)

    ###########
    # Regions #
    ###########
    print("Definition regions...")
    sim.addRegion('CellWall')
    sim.addRegion('Cytoplasm')
    sim.addRegion('Nucleus')
    sim.addRegion('Cajal')
    sim.addRegion('Speckle')
    sim.addRegion('NucEnv')
    sim.addRegion('NPC')
    sim.addRegion('Mito')
    sim.addRegion('Golgi') 
    sim.addRegion('ER') 

    return sim

#########################
# Create the Yeast cell #
#########################
print("Creating Yeast cell...")
cell = lmarray.Cell(initRDME)
cell.setGeometryModel(geometryModel)
cell.setReactionModel(reactionModel)
cell.setDiffusionModel(diffusionModel)
cell.setParticleCounts(particleModel)


# Save simulation
# Geometrical parameters can change, e.g. adding {"nuclSize":[micron(3)]} after file name generate a cell with a nucleus size of 3 um
print("Creating simulation file...")
savedFile = cell.generateLMFiles(filename, {"buildER":[True]})

print("Done. Created:",savedFile)

