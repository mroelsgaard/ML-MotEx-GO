from ase.io import read
import numpy as np
from diffpy.Structure import Structure, Atom
from diffpy.srfit.pdf import PDFContribution, PDFParser, PDFGenerator
from diffpy.srfit.fitbase import FitRecipe, FitResults, Profile, FitContribution
from diffpy.srreal.pdfcalculator import DebyePDFCalculator
from scipy.optimize.minpack import leastsq
import matplotlib.pyplot as plt
def fitting(starting_model, # path for ase.io to read from
            structure_catalogue, # structure catalogue generated with MotEx
            plot, # Bool option to plot or not.
            index,
            qmin=0.7, # Q-min
            qmax=20, # Q-max
            fitrange=(1,10), # r-range for fit
            NumMetals=64, # TODO: structure-gen should be moved out?
            threshold=2.5,
            Experimental_Data="",
           ) -> float:
    """
    This function takes in a 'starting_model', and an 'index' from the 'structure_catalogue'. It generates the 
    corresponding structure and fit it to the 'Experimental_Data'.
    
    returns R-factor from the fit    
    """
    
    # Read structure and divide it into two lists: Atoms we want to iterate (W) and atoms we do not iterate (O)
    stru = read(starting_model)
    xyz = stru.get_positions()
    xyz_Ga = xyz[:NumMetals].copy()
    xyz_O = xyz[NumMetals:len(xyz)].copy()
    keep_O = np.zeros(len(xyz_O))
    h = 0
    
    # Cycle through W atoms and delete W according to index 0's from permutation
    permutations = np.asarray(structure_catalogue)[:,1:]
    for j in range(len(xyz_Ga)):
        if permutations[index][j] == 0:
            xyz_Ga = np.delete(xyz_Ga,j - h,0)
            h = h+1   
            
    # Cycle through all atoms that is not iteratable and test if it is within the threshold distance. Delete atoms with no bonds
    for j in range(len(xyz_O)):        
        for k in range(len(xyz_Ga)):
            dist = np.linalg.norm(xyz_Ga[k] - xyz_O[j])
            if dist < threshold:    
                keep_O[j] = 1
                break
    h = 0            
    for j in range(len(xyz_O)):
        if keep_O[j] == 0:
            xyz_O = np.delete(xyz_O,j - h, 0)
            h += 1
            
    # Create structure for iterable (W) and non-iterable (O) atoms and combine them
    Ga_cluster = Structure([Atom('Ga', xi) for xi in xyz_Ga])
    O_cluster = Structure([Atom('O', xi) for xi in xyz_O])
    cluster = Ga_cluster + O_cluster
    
    # Make a standard cluster refinement using Diffpy-CMI
    # Import the data and make it a PDFprofile. Define the range of the data that will be used in the fit.
    pdfprofile = Profile()
    pdfparser = PDFParser()
    pdfparser.parseFile(Experimental_Data)
    pdfprofile.loadParsedData(pdfparser)
    pdfprofile.setCalculationRange(xmin = fitrange[0], xmax = fitrange[1])

    # Setup the PDFgenerator that calculates the PDF from the structure
    pdfgenerator_cluster = PDFGenerator("G")
    # Add the profile and both generators to the PDFcontribution
    pdfcontribution = FitContribution("pdf")
    pdfcontribution.setProfile(pdfprofile, xname="r") 
    pdfcontribution.addProfileGenerator(pdfgenerator_cluster)
    
    pdfgenerator_cluster.setQmin(qmin)
    pdfgenerator_cluster.setQmax(qmax)
    pdfgenerator_cluster._calc.evaluatortype = 'OPTIMIZED'
    pdfgenerator_cluster.setStructure(cluster, periodic = False)

    # Use scaling factors proportional to molar content
    pdfcontribution.setEquation('mc*G')

    # Define the recipe to do the fit and add it to the PDFcontribution
    recipe = FitRecipe()
    recipe.addContribution(pdfcontribution)

    # Minimize output during fitting 
    recipe.clearFitHooks()

    # Add the scale factor.
    recipe.addVar(pdfcontribution.mc, 1.0, tag = "scale")
    
    # Add the instrumental parameters to the two generators
    pdfgenerator_cluster.qdamp.value = 0.05
    pdfgenerator_cluster.qbroad.value = 0.01

    # Add the delta2 parameters, and make sure it cannot take unphysical values
    recipe.addVar(pdfgenerator_cluster.delta2, 0, name = "delta2_cluster", tag = "delta2")

    # Add ADP and "cell" for the cluster
    phase_cluster = pdfgenerator_cluster.phase
    atoms = phase_cluster.getScatterers()
    lat = phase_cluster.getLattice()

    # lattice "zoomscale" parameters
    recipe.newVar("zoomscale1", 1.0, tag = "lat")
    recipe.newVar("zoomscale2", 1.0, tag = "lat")
    recipe.newVar("zoomscale3", 1.0, tag = "lat")
    recipe.constrain(lat.a, 'zoomscale1')
    recipe.constrain(lat.b, 'zoomscale2')
    recipe.constrain(lat.c, 'zoomscale3')

    # add atomic displacement parameters (adp)
    Ga_cluster = recipe.newVar("Ga_Biso_cluster1", 0.4, tag = 'adp_ga')
    O_cluster = recipe.newVar("O_Biso_cluster1", 0.4, tag = 'adp_o')

    # iterate atomic displacement parameters and constrain
    for atom in atoms:
        if atom.element.title() == "Ga":
            recipe.constrain(atom.Biso, Ga_cluster)
        elif atom.element.title() == "O":
            recipe.constrain(atom.Biso, O_cluster)

    recipe.restrain("zoomscale1", lb = 0.99, ub = 1.01, sig = 0.001)
    recipe.restrain("zoomscale2", lb = 0.99, ub = 1.01, sig = 0.001)
    recipe.restrain("zoomscale3", lb = 0.99, ub = 1.01, sig = 0.001)
    
    #free parameters are set
    recipe.fix('all')
    recipe.free("scale", "lat")

    # Turn off printout of iteration number.
    #recipe.clearFitHooks()

    # We can now execute the fit using scipy's least square optimizer.
    leastsq(recipe.residual, recipe.getValues())
    
    # We calculate the goodness-of-fit, Rwp
    g = recipe.pdf.profile.y
    gcalc = recipe.pdf.evaluate()
    rfactor1 = np.sqrt(sum((g - gcalc)**2) / sum((g)**2))
    
    # if plot is True it will also plot the fit
    if plot:
        print ("FIT RESULTS")
        res1 = FitResults(recipe)
        print (res1)

        # Plot the observed and refined PDF.
        # Get the experimental data from the recipe
        r = recipe.pdf.profile.x
        gobs = recipe.pdf.profile.y

        # Get the calculated PDF and compute the difference between the calculated and measured PDF
        gcalc = recipe.pdf.evaluate()
        baseline = 1.1 * gobs.min()
        gdiff = gobs - gcalc

        # Plot!
        
        fig, ax = plt.subplots()
        ax.plot(r, gobs, 'bo', label="G(r) data")
        ax.plot(r, gcalc, 'r-', label="G(r) fit")
        ax.plot(r, gdiff + baseline, 'g-', label="G(r) diff")
        ax.plot(r, np.zeros_like(r) + baseline, 'k:')
        ax.set(xlabel=r"$r (\AA)$",
               ylabel=r"$G (\AA^{-2})$",
              )
        ax.legend()
        ax.tick_params(direction='in')
        
        # a little bit hacked...
        return fig, ax, rfactor1
    
    # return rwp 
    return rfactor1