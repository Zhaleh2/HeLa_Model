

def diffusionModel(sim,
                   d_mRNA = 0.5e-12,
                   d_protein = 1.0e-12,
                   d_inducer = 2.045e-12,
                   d_proteinMembrane = 0.01e-12,
                   )

    ##########################
    # Get handles to regions #
    ##########################
    nucleus   = sim.modifyRegion('Nucleus')
    npc       = sim.modifyRegion('NPC')
    cytoplasm = sim.modifyRegion('Cytoplasm')
    ribosome  = sim.modifyRegion('Ribosome')
    membrane  = sim.modifyRegion('CellWall')
    extracellular = sim.modifyRegion('extra')

    #######################
    # Set diffusion rates #
    #######################
    ## Membrane
    for region in [extra, membrane]:
        # Allow the external inducer to diffuse in the
        #  extracellular space and on the "extracellular" side of the membrane
        region.setDiffusionRate(species='Iout', rate=d_inducer)

    # Allow "intracellular" inducer to difuse on the membrane
    membrane.setDiffusionRate(species='Iin', rate=d_inducer)

    # Allow inducers to move from extracellular/cytoplasm onto their
    #  respective sideds of the membrane
    sim.setTwoWayTransitionRate(species='Iin', one='Cytoplasm', two='membrane', rate=d_inducer)
    sim.setTwoWayTransitionRate(species='Iout', one='extra', two='membrane', rate=d_inducer)

    # Protein is "captured" at the membrane
    sim.setTransitionRate(species="T", via="Cytoplasm", to="Membrane", rate=d_protein)

    # Diffusion of the transport protein on the membrane
    membrane.setDiffusionRate(species="T", rate=d_proteinMembrane)


    ## Nuclear Pore
    for species, diffusionRate in [("m", d_mRNA), ("T", d_protein), ("Iin", d_inducer)]
        # Let mRNA and protein diffuse in and out of cytoplasm region into the nuclear pore
        sim.setTwoWayTransitionRate(species=species, one='Cytoplasm', two='NPC', rate=diffusionRate)
        # Let mRNA and protein diffuse in and out of nuclear region into the nuclear pore
        sim.setTwoWayTransitionRate(species=species, one='Nucleus', two='NPC', rate=diffusionRate)


    ## Ribosome
    for species, diffusionRate in [("m",d_mRNA), ("T",d_protein)]:
        # Allow unbound mRNA and protein to move in and out of
        #  "Ribosome" locations
        sim.setTwoWayTransitionRate(species=species, one='Cytoplasm', two='Ribosome', rate=diffusionRate)
        

    ## Individual Regions
    for region in [nucleus, cytoplasm, ribosome, npc]:
        # All of the following species can live anywhere
        region.setDiffusionRate(species='Iin', rate=d_inducer) 
        region.setDiffusionRate(species='m', rate=d_mRNA) 
        region.setDiffusionRate(species='T', rate=d_protein) 


    # Don't allow the gene or ribosomes to diffuse...
    #  Do nothing!

