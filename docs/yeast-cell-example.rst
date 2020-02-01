
Example: Yeast Cell Geometry
============================

As a demonstration of how the framework can be used to simulate other
models, the following example will show how a new cell geometry 
can be implemented within the code framework. The geometry 
model is roughly derived from the paper
`Earnest et al. Challenges of Integrating Stochastic Dynamics and Cryo-Electron Tomograms in Whole-Cell Simulations, J. Phys. Chem. B. 2017, 121: 3871-3881 <https://pubs.acs.org/doi/pdf/10.1021/acs.jpcb.7b00672>`_
and roughly approximates the shape of a *S. cerevisiae* cell.

The reaction model for spliceosome assembly presented in the 
manuscript is retained, as this example is primarily used to
show how a new geometrical model can be created. While several 
of the underlying assumptions of the spliceosome model do not
hold in yeast---for instance, the apparent lack of nuclear 
speckles in yeast---the example is illustrative of how a 
new geometry can be simply implemented.

Code for the model can be found in the repository at 
``docs/examples/yeast-cell``.

Updating Simulation Domain
--------------------------
Since the yeast cell from the Earnest et al. manuscript was determined
to be 5.25um x 3.5um with a 1.5 aspect ration, it is much smaller than the
HeLa cell. We change the size of the simulation domain so as not to waste computation 
on blank space. Additionally, the lattice spacing is refined from 64nm
to 28.7nm so that the cell/nuclear walls can be accurately modelled.
These code changes are made in ``yeast.py`` in the ``initRDME`` function:

.. code-block:: python

    def initRDME():
        latticeSpacing = 28.7 # nm
        Dfastest = 6.1e-13 # The fastest diffusion coefficint is required for the calculation of the timestep
        Tsim = 15 # Simulation time in seconds
        sim = RDMESimulation(dimensions=micron(5.51, 5.51, 5.51), # Dimensions correspond to a 5.25um x 3.5um x 3.5um ellipsoidal cell
                             spacing=nm(latticeSpacing), 
                             defaultRegion='extra')
        

Update Geometry Model
---------------------
The first step is to define an ellipsoidal cell shape,
as determiend in the Earnest et al. manuscript. To begin
we define a set of ellipsoidal objects for the cell and
nucleus, respectively along with modified size and quantity
variables for the Cajal bodies, nuclear speckles,  mitochondria
and nuclear pore comples. These numbers were taken directly 
from Earnest et al. or computed as based on the ratio of
yeast:HeLa cell volumes.

.. code-block:: python

    def geometryModel(sim,
                      cellCenter = lm.point(*micron(2.755,2.755,2.755)),
                      cellRadius1 = micron(2.625), # Long dimension of ellipsoid for cell
                      cellRadius2 = micron(1.75), # Short dimension of ellipsoid for cell
                      membraneThickness= nm(28.7), # As determined from the Earnest et al. manuscript
                      nuclSize1 = micron(1.09), # Long dimension of ellipsoid for nucleus
                      nuclSize2 = micron(0.911), # Short dimension of ellipsoid for nucleus
                      speckleRadius=micron(0.05), # Size reduced for smaller nuclear volume
                      cajalRadius=micron(0.07), # Size reduced for smaller nuclear volume
                      poreRadius=micron(0.040), # As determined in Earnest et al.
                      n_cajals=4,
                      n_speckles=20,
                      n_NPCs=139, # As determined in Earnest et al.
                      n_mito=95, # Scaled by ratio of yeast:HeLa cell volume
        ):
    
        ##########################
        # Create cell components #
        ##########################
        # Create an ellipsoid aligned along the X-axis with an
        #  aspect ratio of 1.5 as determined in the Earnest et al. manuscript
        CellWall  = lm.Ellipse(cellCenter, cellRadius1, cellRadius2, cellRadius2, sim.siteTypes['CellWall'])
        Cytoplasm = lm.Ellipse(cellCenter, cellRadius1-membraneThickness, cellRadius2-membraneThickness, cellRadius2-membraneThickness, sim.siteTypes['Cytoplasm'])
    
        # Create an ellipsoid aligned along the X-axis with an
        #  aspect ratio of 1.2 as determined in the Earnest et al. manuscript
        NucEnv    = lm.Ellipsoid(cellCenter, nuclSize1, nuclSize2, nuclSize2 sim.siteTypes['NucEnv'])
        Nucleus   = lm.Ellipsoide(cellCenter, nuclSize1-membraneThickness, nuclSize2-membraneThickness, nuclSize2-membraneThickness, sim.siteTypes['Nucleus'])


Next, we update the Cajal, Speckle and NPC production code to
account for the ellipsoidal nucleus:

.. code-block:: python

    ###### Cajal body: "n_cajal" with radius "cajalRadius" #############################
    dist = [] ; listdist = []
    Cajalpre = []
    accepted = []
    np_accepted = []
    cc = 1 ; point = 0
    while cc < n_cajals + 1:
        x,y,z = np.random.uniform(-0.5,0.5,3)
        Rp = nuclSize2 - cajalRadius - micron(0.1)
        # Account for elliposoidal nature of nucleus in X-axis
        xpp = 1.2*x*Rp + cellCenter.x
        ...

and

.. code-block:: python

    ################ Nuclear Speckles: "n_speckles" with radius "speckleRadius" ##########################
    counter = 0
    taken = 0 
    listdist = [] ; np_acceptspeck = [] ; acceptspeck = [] 
    while counter < n_speckles + 1:
        u = np.random.uniform(-1, 1)
        thet = np.random.uniform(0, 2.0*math.pi)
        rad = np.random.uniform(speckleRadius*1e6 , nuclSize2*1e6 - (speckleRadius*1e6 + 0.1))
        # Account for elliposoidal nature of nucleus
        xm = 1.2*micron(rad*math.sqrt(1-u**2)*np.cos(thet)) + cellCenter.x
        ...

and

.. code-block:: python

    ################# NPC #########################
    NPCpre = []
    taken = 0

    ldist = [] ; npcdist = 0 ; np_npcdist = []
    acceptnpc = []
    npc = 0

    while npc < n_NPCs + 1:
        x,y,z = np.random.uniform(-0.5,0.5,3)
        r = math.sqrt(x**2 + y**2 + z**2)
        # Account for elliposoidal nature of nucleus
        Rp = (nuclSize2 - sim.latticeSpacing)/r
        xpp = 1.2*x*Rp + cellCenter.x
        ...

And then we modify the mitochondria generation code
to account for the smaller cell size:

.. code-block:: python

    ######################## Mitochondria ###############################3
    mit = []
    pp1 = []
    pp2 = []
    Mito = lm.UnionSet(sim.siteTypes['Mito'])
    Mito.thisown=0
    for i in range(1, n_mito+1):
        u = np.random.uniform(-1, 1)
        thet = np.random.uniform(0, 2.0*math.pi)
        rad = np.random.uniform(0.5, 2.5) # Modified for smaller cell
        xm = 1.5*micron(rad*math.sqrt(1-u**2)*np.cos(thet)) + micron(9) 
        ...

We remove the definition of the Golgi apparatus, as it was not
modelled in Earnest et al.

.. code-block:: python

    # Add all geometries to the simulation
    sim.lm_builder.addRegion(CellWall)
    sim.lm_builder.addRegion(Cytoplasm)
    sim.lm_builder.addRegion(NucEnv)
    sim.lm_builder.addRegion(Nucleus)
    sim.lm_builder.addRegion(NPC)
    sim.lm_builder.addRegion(Cajal)
    #sim.lm_builder.addRegion(Cajal) # Removed
    sim.lm_builder.addRegion(Speckle)
    sim.lm_builder.addRegion(Mito)

Finally we remove the ER function ``createERCellularAutomaton`` from ``yeast_geometry.py``
as the ER is ignored in the manuscript.

Final Notes
-----------
While the above code is illustrative of how the a different 
geometry could be modeled within the framework presented here, we note
that the reaction model is taken wholesale from our manuscript.
As mentioned above, several of the assumptions, such as presence
of nuclear speckles and the separation of spliceosome biogenesis
across the nuclear membrane, may not be applicable in yeast.
Regardless, we believe this code is illustrative of how a 
different geometry model could be easily implemented within the framework.

