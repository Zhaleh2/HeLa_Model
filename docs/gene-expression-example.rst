
Example 1: Inducible Genetic Swith
==================================

As a demonstration of how the framework can be used to simulate other
models, the following example will show how another reaction model
can be implemented within the HeLa cell is shown. The reaction
model is derived from the paper
`Earnest et al. Challenges of Integrating Stochastic Dynamics and Cryo-Electron Tomograms in Whole-Cell Simulations, J. Phys. Chem. B. 2017, 121: 3871-3881 <https://pubs.acs.org/doi/pdf/10.1021/acs.jpcb.7b00672>`_
and roughly approximates the galactose inducible genetic circuit 
from *S. cerevisiae*.

The model consists of a switchable gene, transcription of the 
coded mRNA, and translation of a transport protein at the ribosome.
To allow for switching behavior, the reaction system contains
transport terms for the galactose molecules which induce the 
switch from an "off" state to an "on" state.

In the following code snippets, we demonstrate how the HeLa
cell code can be easily modified to model this system within the
HeLa cell geometry presented in the manuscript. Code for the model
can be found in the repository in ``docs/examples/gene-expression-example``.

Specifying the Species
----------------------
In specifying a new model, the first modifications that need to
be made is the specification of the driver file ``hela.py``. As shown
in the following code snippet, the species that make up the system are 
defined.

.. code-block:: python

    ###################
    # Specify Species #
    ###################
    print("Defining Species...")
    # Iout - External inducer
    # Iin - Internal inducer
    # Doff - DNA in the non-transcribing state
    # Don - DNA in a transcribing state
    # m - mRNA
    # R - Ribosome
    # R:m - Ribosome bound to mRNA
    # T - Transporter
    # T:I - Transporter bound to inducer
    species = ["Iout","Iin", "Doff", "Don", "m", "R", "R:m", "m'", "T", "T:I"]


One addition needs to be made to the spatial model to support the
reaction model as described in the manuscript; namely, the addition
of specific locations that represent fixed ribosomes:

.. code-block:: python

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
    sim.addRegion('Ribosome') # Added

The above addition represents a new "site type" that can be used
to confine reactions to a particular spatial location within
the simulated cell. Because ribosomes are immobile in the manuscript,
we will use this new region as fixed location to place ribosome 
"particles", which can catalze the translation reactions.


Update the Geometry Model to Add "Ribosome" Sites
-------------------------------------------------
After specifying the spatial model, we can add ribosome sites to
the cell geometry. First, a new parameter is added to the 
geometry model function that allows for the number of ribosomes
to be easily modified. The function definition is modified in 
``hela_geometry.py`` as:

.. code-block:: python

    def geometryModel(sim,
                      cellCenter = lm.point(*micron(9,9,9)),
                      cellRadius = micron(8.9),
                      membraneThickness= micron(0.128),
                      nuclSize= micron(4.15),
                      speckleRadius=micron(0.35),
                      cajalRadius=micron(0.5),
                      poreRadius=micron(0.083),
                      n_cajals=4,
                      n_speckles=20,
                      n_NPCs=1515,
                      n_mito=2000,
                      n_ribo=750, # Added Line
                      fb = 0.9,
                      fd = 0.8,
                      steps=11,
                      lsize=288,
                      buildER=True,
                      limits=[lambda x: x<= 65.0**2,lambda x: x> 139.0**2]
        ):


Next, code is added at the end of the function to randomly
place ribosome locations within the cytoplasm:

.. code-block:: python

    ######################### Ribosomes #########################################
    # Add 750 ribosomes to the cell's cytoplasm randomly
    if n_ribo > 0:
        added = 0
        for x,y,z in np.random(5*n_ribo*3).reshape((5*n_ribo,3))*cellRadius:
            if Cytoplasm.intersects(lm.Sphere(lm.Point(x,y,z)), nm(32), sim.siteTypes['Cytoplasm']):
                sim.setLatticeSite((x,y,z), 'Ribosome')
                added += 1
                
            if added >= n_ribo:
                break


This code could be modified to allow explicit locations be 
specified if additional information is available, perhaps
from cryo-electron tomography as was done in the manuscript.


Update Reaction Model
---------------------
Once the spatial and geometry models are specified, the reaction
model can be specified. Here the reaction model of the HeLa
cell is completely replaced with the by that presented in the
paper.

To allow simple tuning of the model, rate constants are specified
as parameters to the ``reactionModel`` function. Subsequently,
reactions are speficied in the spatial regions where they take
place.


.. code-block:: python

    def reactionModel(sim,
                  kgnOn=1.599,
                  kts=6.202e-3,
                  ktlInit=7.043e-3,
                  ktlTerm=1.393,
                  kmDcy=7.889e-4,
                  kts_prime=5.895e-5,
                  ktlTerm_prime=1.101,
                  kmDcy_prime=5.776e-4,
                  ktxDif=2.33e-3,
                  ktxOn=2.134,
                  ktx=12.0,
                  ktxOff=0.12,
                  ktDcy=2.567e-4,
                  ):


        ##########################
        # Adjust rates by volume #
        ##########################
        scale = 157863.12  # = 6.022e23 * 64^3 e-24

        ##########################
        # Get handles to regions #
        ##########################
        nucleus   = sim.modifyRegion('Nucleus')
        cytoplasm = sim.modifyRegion('Cytoplasm')
        ribosome  = sim.modifyRegion('Ribosome') # Added
        membrane  = sim.modifyRegion('CellWall') # Added

        #####################
        # Add all Reactions #
        #####################
        nucleus.addReaction(reactant=('Iin','Doff'), product=('Don'), rate=kgnOn) # inducer/TF binding
        nucleus.addReaction(reactant=('Don'), product=('Don','m'), rate=kts) # transciption
        ribosome.addReaction(reactant=('R','m'), product=('R:m'), rate=ktlInit) # SSU/mRNA association
        ribosome.addReaction(reactant=('R:m'), product=('R','m','T'), rate=ktlTerm) # translation elongation
        for region in [nuclus,ribosome]:
            region.addReaction(reactant=('m'), product='', rate=kmDcy) # mRNA degradation
        ribosome.addReaction(reactant=('R:m'), product=('R'), rate=kmDcy) # mrnaDegradation
        nucleus.addReaction(reactant='', product='m', rate=kts_prime) # transciption (other)
        ribosome.addReaction(reactant=('R','m_prime'), product=('R:m_prime'), rate=ktlInit) # SSU/mRNA association (other)
        ribosome.addReaction(reactant=('R:m_prime'), product=('R','m_prime'), rate=ktlTerm) # translation elongation (other)
        for region in [nuclus,ribosome]:
            region.addReaction(reactant=('m_prime'), product='', rate=kmDcy_prime) # mRNA degradation (other)
        ribosome.addReaction(reactant=('R:m_prime'), product=('R'), rate=kmDcy_prime) # mRNA degradation (other)
        membrane.addReaction(reactant=('Iin'), product=('Iout'), rate=ktxDif) # passive diffusional transport
        membrane.addReaction(reactant=('Iout'), product=('Iin'), rate=ktxDif) # passive diffusional transport
        membrane.addReaction(reactant=('T','Iout'), product=('T:I'), rate=ktxOn) # transporter/inducer association
        membrane.addReaction(reactant=('T:I'), product=('T','Iin'), rate=ktx) # active inducer transport
        membrane.addReaction(reactant=('T:I'), product=('T','Iout'), rate=ktxOff) # transporter/inducer dissociation
        for region in [cytoplasm, membrane]:
            region.addReaction(reactant=('T'), product='', rate=ktDcy) # transporter degradation
            region.addReaction(reactant=('T:I'), product='', rate=ktDcy) # transporter degradation
        print("Reactions are set!")


Update Particle Model
---------------------
...

.. code-block:: python

    def particleModel(sim,
                      n_ribo=750,
                      n_Doff=1,
                      n_Don=0,
                      n_m=1,
                      n_T=100,
                      n_I=1000000,
                      ):
    
        ####################
        # Add all Species #
        ###################
    
        sim.addParticles(species='R', region='Ribosome', count= n_ribo)
        sim.addParticles(species='T', region='CellWall', count= n_T)
        sim.addParticles(species='Doff', region='Nucleus', count= n_Doff)
        sim.addParticles(species='Don', region='Nucleus', count= n_Don)
        sim.addParticles(species='m', region='Nucleus', count=n_m)
        sim.addParticles(species='Iout', region='extra', count=n_I)
    
        print("Particles were added!")


Update Diffusion Model
----------------------
TODO

Final Notes
-----------
While the above code is illustrative of how the a different 
reaction model could be modeled within the HeLa cell, we note
that the model is implemented as presented in the paper, which
is tailored for a cell the size of an average yeast microbe.
A proper treament would require modifiying particle numbers
to account for the differences in size between a yeast and a HeLa
cell, and diffusion and reaction rates that are appropriate
for a genetic switch in a Human cell. That said, it exemplifies
how the code can be easily modified to model processes
not originally anticipated in our manuscript.

