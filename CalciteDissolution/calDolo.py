import importlib.util
import os

description = "Test gems, calcite dissolution"
timeUnit = 's'

maxTime = 21e3
timeStep = 500

chemModule = 2

systemSpeciation = ["Ca(CO3)@","Ca(HCO3)+","Ca+2","CaOH+","Mg(CO3)@","Mg(HCO3)+","Mg+2","MgOH+","CO2@","CO3-2","HCO3-","CH4@","ClO4-","Cl-","H2@","O2@","OH-","H+","H2O@","CO2","CH4","H2","O2","Cal","Dis-Dol","Sn"]

nonTrivialDC = {"Cal": {"Ca" : 1, "C" : 1, "O" : 3},
                    "Dis-Dol": {"Mg" : 1, "C" : 2, "Ca" : 1, "O" : 6}
                    }

chemPath ="/volatile/home/al274877/PICCTS/Resources/CalciteIC/CalciteIC-dat.lst"
trsptPath = "calDoloGems.mph"

independentComponents = ["C","Ca","Cl","H","Mg","O","Sn","Zz"] # attention to the order ! should match the order in the .lst files

# initialConditions = "/volatile/home/al274877/PICCTS/CalDolo/ic.txt"

firstStepEquilibrium = True
PIDnbr = 2 # xGEMS only, COMSOL is natively parallelised


enginePath  ="" # change your path here
if __name__ == "__main__":
    spec = importlib.util.spec_from_file_location(rf"{enginePath}/engine",rf"{enginePath}/engine.py")
    with open("store.txt", "w", encoding="utf-8") as fichier:
        fichier.write(os.path.dirname(os.path.abspath(__file__))+'\n')
        fichier.write(os.path.basename(__file__))
    code = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(code)
    code.main()
