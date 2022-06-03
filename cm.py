# routines to output CrystalMaker files
import numpy as np
import matplotlib as mpl
import matplotlib.cm
import random

def Make_CrystalMakerFile(elements, 
                            xyz, 
                            AtomContributionValues, 
                            m, 
                            saveResults, 
                            threshold,
                            num_metals: int):
    """
     Creates a CrystalMaker file with color-coded atoms based on the contribution
     values that are output from ML MotEx

     Inputs:
      elements
      AtomContributionValues
      m
      saveResults
      threshold
      num_metals
    """


    # Output a crystalmaker file header to visualize the results
    CrystalMaker = open(saveResults+'CrystalMaker_MotEx.cmtx', 'w')

    CrystalMaker.write("MOLE  CrystalMaker molecule format\n")
    CrystalMaker.write("TITL  Molecule\n\n")
    CrystalMaker.write("! Model type\n")
    CrystalMaker.write("MODL  1\n\n")

    CrystalMaker.write("! Depth fading settings\n")
    CrystalMaker.write("DCUE  1.000000 0.212899 0.704686\n\n")

    CrystalMaker.write("! Colour definitions:\n")
    CrystalMaker.write("TYPE\n")

    # Assign colors to all the atoms
    for iter, element in enumerate(elements):
        if iter < num_metals:
            CrystalMaker.write(element + str(iter+1) + " 1.32 ")
            rgb1 = m.to_rgba(AtomContributionValues[iter])[:-1][0]
            rgb2 = m.to_rgba(AtomContributionValues[iter])[:-1][1]
            rgb3 = m.to_rgba(AtomContributionValues[iter])[:-1][2]
            CrystalMaker.write(f'{rgb1} {rgb2} {rgb3}')
            CrystalMaker.write("\n")
        else:
            CrystalMaker.write(f'{element}{iter+1} 0.66')
            rgb1 = mpl.colors.to_rgb("#FF0000")[0]
            rgb2 = mpl.colors.to_rgb("#FF0000")[1]
            rgb3 = mpl.colors.to_rgb("#FF0000")[2]
            CrystalMaker.write(f'{rgb1} {rgb2} {rgb3}')
            CrystalMaker.write("\n")
    
    CrystalMaker.write("\n")
    CrystalMaker.write("! Atoms list\n")
    CrystalMaker.write("! Bond Specifications\n")
    
    # Assign bonds between the atoms
    for iter, element in enumerate(elements):
        if iter < num_metals:
            NI_elements = np.delete(np.unique(elements), np.where(np.unique(elements) == element)[0])
            for NI_element in NI_elements:
                CrystalMaker.write(f"BMAX {element}} {NI_element}  {threshold}"))
                CrystalMaker.write("\n")
    
    CrystalMaker.write("\n")
    CrystalMaker.write("! Atoms list\n")
    CrystalMaker.write("ATOM\n")
    
    # Assign coordinates to the atoms
    for iter, element in enumerate(elements):
        if iter < num_metals:
            CrystalMaker.write(f'{element} {element} {iter+1} {xyz[iter][0]} {xyz[iter][1]} {xyz[iter][2]}\n')
        else:
            CrystalMaker.write(f'{element} {element} {iter+1} {xyz[iter][0]} {xyz[iter][1]} {xyz[iter][2]}\n')

    CrystalMaker.close()
    
    return None
