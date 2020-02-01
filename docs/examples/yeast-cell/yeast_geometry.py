import lm
from pyLM.units import *
import numpy as np
import scipy as sp
from scipy import ndimage
import math
import itertools


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


    # The following code snippets are required to prevent python from 
    #  garbage collecting the objects when the function exits, allowing
    #  lattice microbes to refer to them later when building the cell
    CellWall.thisown = 0
    Cytoplasm.thisown = 0
    NucEnv.thisown = 0
    Nucleus.thisown = 0

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
        ypp = y*Rp + cellCenter.y  
        zpp = z*Rp + cellCenter.z 
        caj = [xpp,ypp,zpp]
        listdist = []
        point += 1
        if (point == 1 ):
            accepted.append(caj)
            np_accepted = np.array(accepted)
        for i in range (0, cc):
            dist =  np.sqrt(np.sum((np.array(caj)-np_accepted[i])**2))
            listdist.append(dist) 
        np_dist = np.array(listdist)
        if(np.all(np_dist > 2*cajalRadius + micron(0.05))):
            accepted.append(caj)
            np_accepted = np.array(accepted)
            listdist = []
            cc += 1
    NCajal = np.array(accepted)

    cc = [] 
    Cajal = lm.UnionSet(sim.siteTypes['Cajal'])
    Cajal.thisown = 0
    for k in range(0, n_cajals):
        w = lm.Sphere(lm.point(NCajal[k][0],
                               NCajal[k][1],
                               NCajal[k][2]),
                      cajalRadius,
                      sim.siteTypes['Cajal'])
        w.thisown=0
        cc.append(w)
        Cajal.addShape(w)
    print("Cajal done!")

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
        ym = micron(rad*math.sqrt(1-u**2)*np.sin(thet)) + cellCenter.y
        zm = micron(rad*u) + cellCenter.z
        point = [xm, ym, zm]
        listdist = [] ; dist = 0 ; np_dist = []
        taken += 1
        if ( taken == 1 ):
            acceptspeck.append(point)
            np_acceptspeck = np.array(acceptspeck)
        for i in range (0, counter):
            dist =  np.sqrt(np.sum((point-np_acceptspeck[i])**2))
            listdist.append(dist)
        np_dist = np.array(listdist)
        tt = 0
        if(np.all(np_dist > 2*speckleRadius + micron(0.1))):
            for k in range (0, n_cajals):
                if (np.sqrt(np.sum((np.array(point)-NCajal[k])**2)) > speckleRadius + cajalRadius + micron(0.1)):
                    tt += 1
            if (tt == 4):
                acceptspeck.append(point)
                np_acceptspeck = np.array(acceptspeck)
                counter += 1
                listdist = []
    NSpeck = np.array(acceptspeck)

    collection = []
    j = 0
    Speckle = lm.UnionSet(sim.siteTypes['Speckle'])
    Speckle.thisown = 0
    for k in range(0, n_speckles):
        q = lm.Sphere(lm.point(NSpeck[k][j],
                               NSpeck[k][j+1],
                               NSpeck[k][j+2]),
                      speckleRadius,
                      sim.siteTypes['Speckle'])
        q.thisown=0
        collection.append(q)
        Speckle.addShape(q)
    print("Speckles done!")

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
        ypp = y*Rp + cellCenter.y
        zpp = z*Rp + cellCenter.z
        pores = [xpp,ypp,zpp]
        NPCpre.append(pores)
        ldist = [] ; np_npcdist = []
        taken += 1
        if (taken == 1 ):
            acceptnpc.append(pores)
            np_acceptnpc = np.array(acceptnpc)
        for i in range (0, npc):
            npcdist =  np.sqrt(np.sum((np.array(pores)-np_acceptnpc[i])**2))
            ldist.append(npcdist)
        np_npcdist = np.array(ldist)
        if(np.all(np_npcdist > 0.2e-6)):
            acceptnpc.append(pores)
            np_acceptnpc = np.array(acceptnpc)
            ldist = []
            npc += 1

    NP = np.array(NPCpre)

    collect = []
    j = 0
    NPC = lm.UnionSet(sim.siteTypes['NPC'])
    NPC.thisown = 0
    for k in range(0, n_NPCs):
        s = lm.Sphere(lm.point(NP[k][j],
                               NP[k][j+1],
                               NP[k][j+2]),
                      poreRadius,
                      sim.siteTypes['NPC'])
        s.thisown=0
        collect.append(s)
        NPC.addShape(s)

    print("NPCs done!")

    ####################### Golgi Apparatus  #################################################
    # Skip golgi apparatus as done in Earnest et al.

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
        ym = micron(rad*math.sqrt(1-u**2)*np.sin(thet)) + micron(9) 
        zm = micron(rad*u) + micron(9)
        centers = [xm, ym, zm]
        if i < 666 :
            point1 = [xm - micron(0.125), ym-micron(0.125) , zm]
            point2 = [xm + micron(0.125), ym+micron(0.125) , zm]
        elif  i >= 666 and i < 1332 :
            point1 = [xm, ym - micron(0.125), zm- micron(0.125)]
            point2 = [xm, ym + micron(0.125), zm+ micron(0.125)]
        elif  i >= 1332 :
            point1 = [xm - micron(0.125), ym , zm - micron(0.125)]
            point2 = [xm + micron(0.125), ym , zm + micron(0.125)]
            
        pp1.append(point1)
        pp2.append(point2)
    mp1 = np.array(pp1)
    mp2 = np.array(pp2)

    for h in range(0, n_mito):
        m = lm.Capsule(lm.point(mp1[h][0], mp1[h][1], mp1[h][2]),
                       lm.point(mp2[h][0], mp2[h][1], mp2[h][2]),
                       micron(0.249),
                       sim.siteTypes['Mito'])
        m.thisown = 0
        mit.append(m)
        Mito.addShape(m)
    print("Mito done!")



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

    ######################### ER #########################################
    # Ignored in Earnest et al. as it occupies 2.2% of cell volume
 
