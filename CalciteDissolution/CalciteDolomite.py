import importlib.util
import os

description = "Calcite dissolution and dolomite precipiation. Test case presented by Kulik et al. 2013, Azad et al. 2016, ..."
timeUnit = 's'

maxTime = 21e3
timeStep = 500
PIDnbr = 1 # number of PID(s) carrying xGEMS calculation
firstStepEquilibrium = True # this keyword will first equilibrate your initial conditions.

chemModule = 'xgems'

systemSpeciation = ["Ca(CO3)@","Ca(HCO3)+","Ca+2","CaOH+","Mg(CO3)@","Mg(HCO3)+","Mg+2","MgOH+","CO2@","CO3-2","HCO3-","CH4@","ClO4-","Cl-","H2@","O2@","OH-","H+","H2O@","CO2","CH4","H2","O2","Cal","Dis-Dol","Sn"]

nonTrivialDC = {"Cal": {"Ca" : 1, "C" : 1, "O" : 3},
                    "Dis-Dol": {"Mg" : 1, "C" : 2, "Ca" : 1, "O" : 6}
                    }

chemPath ="Resources/CalciteIC/CalciteIC-dat.lst" # change this path if necessary
trsptPath = "CalciteDolomite.mph"

independentComponents = ["C","Ca","Cl","H","Mg","O","Sn","Zz"] # attention to the order ! should match the order in the .lst files

initialConditions = "CalciteDissolutionIC.txt"

firstStepEquilibrium = True

enginePath = "../Source" # change your path here
if __name__ == "__main__":
    spec = importlib.util.spec_from_file_location(rf"{enginePath}/engine",rf"{enginePath}/engine.py")
    with open("store.txt", "w", encoding="utf-8") as fichier:
        fichier.write(os.path.dirname(os.path.abspath(__file__))+'\n')
        fichier.write(os.path.basename(__file__))
    code = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(code)
    code.main()
