import importlib.util
import os

maxTime = 72000
timeStep = 720

chemModule = 'orchestra' # you can change with 'phreeqc'.
operatorSplitting = 'snia' # change the operator splitting using 'additive', 'alternative', 'strang', 'symmetrical' and 'snia'
firstStepEquilibrium = True # this keyword will first equilibrate your initial conditions.

# communication framework is identical between PhreeqC and ORCHESTRA.
# So PhreeqC and ORCHESTRA share the same speciation, database, initial conditions and COMSOL file.
# Only two keywords truly depend on the module chosen (speciesAttributes and supplementarySolution), due to the different formalism of the solvers
systemSpeciation = ["CaX2", "KX", "NaX", "NH4X", "Ca+2", "CaOH+","Cl-", "H+", "H2", "K+", "N2", "Na+", "NaOH", "NH3", "NH4+", "NO2-", "NO3-", "O2", "OH-"]
trsptPath = "CationExchange.mph"
primarySpeciesAq = ["Ca","Cl", "K", "N", "Na"] # Will soon be retrieved from the database ...


if chemModule == 'orchestra':
    initialConditions = "orchestraIC.txt"
    chemPath = 'chemistry1.inp' # Note that 'object2025.txt' must be located in the same directory as this file
    
    speciesAttributes = {
        "con" : [spc for spc in systemSpeciation if spc not in ["CaX2", "KX", "NaX", "NH4X"]],
        "solid" : {'CaX2':'Ca','KX':'K','NaX':'Na','NH4X':'NH4'},
        }
else:
    initialConditions = "phreeqcIC.txt"
    chemPath = "phreeqc.dat" 

    speciationEch = ["CaX2", "KX", "NaX", "NH4X"] # Will soon be retrieved from the database ...
    
    supplementarySolution = 'pe 12.5 O2(g) -0.68'

enginePath = r"...\Source"
if __name__ == "__main__":
    spec = importlib.util.spec_from_file_location("engine", os.path.join(enginePath, "engine.py"))
    with open("store.txt", "w", encoding="utf-8") as fichier:
        fichier.write(os.path.dirname(os.path.abspath(__file__))+'\n')
        fichier.write(os.path.basename(__file__))
    code = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(code)
    code.main()