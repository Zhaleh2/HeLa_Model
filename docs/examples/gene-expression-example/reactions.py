

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
    npc       = sim.modifyRegion('NPC')
    speckle   = sim.modifyRegion('Speckle')
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
