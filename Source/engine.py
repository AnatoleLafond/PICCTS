import time
PICCTSstart = time.time()
import pandas as pd
import sys 
import importlib.util
import os
import shutil

import re

from pathlib import Path
current_dir = Path(__file__).parent
if str(current_dir) not in sys.path:
    sys.path.insert(0, str(current_dir))

with open("store.txt", 'r', encoding='utf-8') as fichier:
    inputPath =  fichier.readline().strip() 
    nameInput =  fichier.readline().strip()


module_name = os.path.splitext(nameInput)[0]
file_path = os.path.join(inputPath, nameInput)
spec = importlib.util.spec_from_file_location(module_name, file_path)

PICCTS_input = importlib.util.module_from_spec(spec)
spec.loader.exec_module(PICCTS_input)

maxTime = getattr(PICCTS_input, 'maxTime', 1)
timeStep = getattr(PICCTS_input, 'timeStep', 1)
stepReprise = getattr(PICCTS_input, 'stepReprise', 0)
dtPICCTS = getattr(PICCTS_input, 'dtPICCTS', [])
if dtPICCTS and dtPICCTS[0]==0: del dtPICCTS[0]

if not dtPICCTS:
    dtPICCTS = [timeStep]
    while dtPICCTS[-1] < maxTime:
        dtPICCTS += [dtPICCTS[-1] + timeStep]

if stepReprise:
    # stepReprise was the last completed step. The Xeme step is the (X-1)eme step for piccts
    # So the Xeme step is the step to be completed
    startTimeStep = start = stepReprise
    timeStepReprise = dtPICCTS[stepReprise-1]
else:
    timeStepReprise = 0
    startTimeStep = start = 0


comsolData = getattr(PICCTS_input, 'comsolData', ['data1'])
outputCopy = getattr(PICCTS_input, 'outputCopy', True)
sortiesDossier= ["primSpecies","outputSpeciation","outputHistory"]
for data in comsolData:
    sortiesDossier += [f'outputCOMSOL{data}',f'VTU{data}']

paths = {} # variables are keys and paths are values
if getattr(PICCTS_input, 'intermediateOutput', True) :
    for i,doss in enumerate(sortiesDossier):
            paths[doss] = os.path.join(inputPath, sortiesDossier[i])
            if os.path.exists(paths[doss]) and not stepReprise:
                shutil.rmtree(paths[doss]) # Quite dangerous but I like risk (will permanently delete your previous files when running a new PICCTS run)
                os.makedirs(paths[doss])
            elif not stepReprise:
                os.makedirs(paths[doss])

    outputName = ["trsptPath","chemPath","initialConditions"]
    outputDefault = ["comsol.mph","phreeqc.dat","ic.txt"]
    for i, txt in enumerate(outputName):
        paths[txt] = getattr(PICCTS_input, txt, None)
        if paths[txt] == None :
            paths[txt] = getattr(PICCTS_input, txt, outputDefault[i])


def main():
    def writeTime(tps, arr=10):
        if tps >= 3600 * 24:
            return f"{tps / (3600 * 24):.{arr}f} d"
        elif tps >= 3600:
            return f"{tps / 3600:.{arr}f} h"
        elif tps >= 60:
            return f"{tps / 60:.{arr}f} min"
        else:
            return f"{tps:.{arr}f} sec"
    
    def readInputFile(txtPath, coord,  inputHeaders = None) :
        dataframe = pd.read_csv(txtPath,sep=r"\s+",comment="%",dtype=float)
        if (coord + inputHeaders) != dataframe.columns.tolist():
            dataframe.columns = coord + inputHeaders
        return dataframe
    
    # 
    
    def extract_master_species(filepath):
        '''
        will extract the primary species from a phreeqc formalism database
        so what the user does not have anymore to report the different primary species
        To continue ...
        '''

        masterSpecies = {
            'SOLUTION_MASTER_SPECIES': [],
            'SURFACE_MASTER_SPECIES': [],
            'EXCHANGE_SPECIES': [],
        }

        stop_keywords = {
            "SOLUTION_SPECIES",
            "SURFACE_SPECIES",
            "EXCHANGE_MASTER_SPECIES",
            "EQUILIBRIUM_PHASES",
            "REACTION",
            "KINETICS",
            "RATES",
            "END"
        }

        current_block = None
        with open(filepath, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"): continue
                upper_line = line.upper()
                if upper_line in masterSpecies:
                    current_block = upper_line
                    continue
                if current_block and any(k in upper_line for k in stop_keywords):
                    current_block = None
                    continue
                if current_block:
                    if line.lower().startswith("log_k"): continue
                    if "=" in line:
                        lhs, rhs = line.split("=", 1)
                        lhs_last = lhs.strip().split()[-1]
                        rhs_first = rhs.strip().split()[0]
                        if lhs_last == rhs_first: continue
                        if not rhs_first.replace(".", "").isdigit(): masterSpecies[current_block].append(rhs_first)
                    parts = line.split()
                    if parts and current_block != "EXCHANGE_SPECIES": masterSpecies[current_block].append(parts[0])

        return masterSpecies
    

    
    
        
    # Decompose secundary species into primary species (e.g. Na2S2O3 into 2Na, 2S and 3O)
    def decomposingIntoPrimSpecies(formula, primarySpecies,):
        def multiply_dict(d, factor): 
            return {k: v * factor for k, v in d.items()}
        def merge_dicts(a, b): 
            for k, v in b.items():
                a[k] = a.get(k, 0) + v
            return a
    
        charge_match = re.search(r'([+-]\d*)$', formula)
        if charge_match:
            charge = charge_match.group(1)
            formula = formula[:charge_match.start()]
        else:
            charge = None
    
        formula = re.sub(r'__([0-9]+)', r')\1', formula) 
    
        primaryspecies_trié = sorted(primarySpecies, key=len, reverse=True)
        pattern = re.compile('|'.join(re.escape(ps) for ps in primaryspecies_trié)) 
    
        index = 0
        tokens = []
    
        while index < len(formula): 
            m = pattern.match(formula, index) 
            if m:
                name = m.group(0)
                j = index + len(name) 
                coef_match = re.match(r'\d+', formula[j:]) 
                if coef_match:
                    count = int(coef_match.group())
                    j += len(coef_match.group()) 
                else:
                    count = 1 
                tokens.append(('group', name, count)) 
                index = j 
                continue
    
            if formula[index] == '(':
                tokens.append(('(',))
                index += 1
                continue
            elif formula[index] == ')': 
                j = index + 1
                while j < len(formula) and formula[j].isdigit():
                    j += 1
                multiplicateur = int(formula[index+1:j]) if j > index+1 else 1 
            
                tokens.append((')', multiplicateur))
                index = j
                continue
    
            index += 1
        stack = []
        current = {}
    
        for token in tokens:
            if token[0] == 'group': 
                name, count = token[1], token[2]
                current[name] = current.get(name, 0) + count 
            elif token[0] == 'element':
                elem, count = token[1], token[2]
                current[elem] = current.get(elem, 0) + count
            elif token[0] == '(':
                stack.append(current)
                current = {}
            elif token[0] == ')':
                multiplicateur = token[1]
                current = multiply_dict(current, multiplicateur)
                prev = stack.pop()
                current = merge_dicts(prev, current)
                
        return current, charge
   
    
    completeCpldParam = [
        ['viscos','viscosity(mPa·s)','v','visc'],
        ['tk','temperature(K)','t','temp','tempK','tK','temperatureK','temperature'],
        ['tc','temperature(°C)','tC','tempC','temperatureC'],
        ['mu','ionicStrength(mol/kgw)','IS','is'],
        ['EDL','EDL_','edl'],
        ['TOT','TotAq_','tot_aq'],
        ]
    coupledParameters = {} 
    partialCoupledParameters = getattr(PICCTS_input, 'partialCoupledParameters', [])
    if getattr(PICCTS_input, 'coupledParameters', None):
        for cpld in getattr(PICCTS_input, 'coupledParameters', None): 
            for cpldList in completeCpldParam:
                if isinstance(cpld, dict): 
                    for key in cpld.keys():
                        if key in cpldList:
                            val = cpld[key]
                            if all(isinstance(v, list) for v in val):
                                coupledParameters[cpldList[1]] = [cpldList[0]] + val
                            else:
                                coupledParameters[cpldList[1]] = [cpldList[0], val]
                                break
                elif cpld in cpldList:
                    coupledParameters[cpldList[1]] = [cpldList[0]]
                    break
                elif cpldList == completeCpldParam[-1]:
                    partialCoupledParameters += [cpld]
                
    coupledParametersNames = []
    for name in coupledParameters.keys():
        val = coupledParameters[name]
        if len(val) > 1 and isinstance(val[1], list):
            for spc in val[1]:
                coupledParametersNames += [f"{name}{spc}"]
        else: coupledParametersNames += [f"{name}"]
    maillesChargeGeom = []
    specieChargeGeom = []
    if getattr(PICCTS_input, 'speciesChargeGeometry', None):
        for (maille, spc) in getattr(PICCTS_input, 'speciesChargeGeometry', None): 
            maillesChargeGeom += [maille]
            specieChargeGeom += [spc]
    
    phases = getattr(PICCTS_input, "phases", "")
    if getattr(PICCTS_input, "fixpH", None):
        phases += "Fix_ph\nH+=H+; log_k 0"
    
    centralDict = { # gather keywords which do not depend on components
        "commMtrx": readInputFile(paths['initialConditions'],['x', 'y', 'z'][:getattr(PICCTS_input, 'geometry', 1)],(partialCoupledParameters + getattr(PICCTS_input, 'systemSpeciation', None))), # pd.read_csv(paths['IC'],sep=r"\s+",comment="%",dtype=float),
        "warningRun": 0,
        "AcidicEcho": pd.DataFrame(),
        "waitingTime" : 0,
        "couplingInfo" : [ ['SNIA','Strang','Alternative','Additive','Symmetrical'][getattr(PICCTS_input, 'operatorSplitting', 1)-1],
                          ['PhreeqC','xGEMS','ORCHESTRA'][getattr(PICCTS_input, 'chemModule', 1)-1],
                          ['COMSOL'][getattr(PICCTS_input, 'trsptModule', 1)-1],
            ],
        "paths" : paths,
        "system" : getattr(PICCTS_input, 'system', 1),
        "stepReprise" : getattr(PICCTS_input, 'stepReprise', False),
        "PIDnbr": getattr(PICCTS_input, 'PIDnbr', 1),
        "rnvllmt": getattr(PICCTS_input, 'renouvellement', None),
        "systemSpeciation" : getattr(PICCTS_input, 'systemSpeciation', None),
        "timeUnit" : getattr(PICCTS_input, 'timeUnit', 's'),
        "geometry" : getattr(PICCTS_input, 'geometry', 1),
        "coupledParameters" :  coupledParameters,
        "coupledParametersNames" : coupledParametersNames,
        "partialCoupledParameters" : partialCoupledParameters,
        "coord" : ['x', 'y', 'z'][:getattr(PICCTS_input, 'geometry', 1)],
        "chemModule" : getattr(PICCTS_input, 'chemModule', 1),
        "trsptModule" : getattr(PICCTS_input, 'trsptModule', 1),
        "firstStepEquilibrium" : getattr(PICCTS_input, 'firstStepEquilibrium', False),
        }
    
    if centralDict["couplingInfo"][1] == 'PhreeqC':
        import phreeqc
        speciationLauncher = {
            1: phreeqc.spct,
        }
        
        centralDict.update({
        "primToSecSpecies" : getattr(PICCTS_input, 'primToSecSpecies', {}),

        "chemPath" : getattr(PICCTS_input, 'chemPath', 'phreeqc.dat'),
        "current" : getattr(PICCTS_input, 'current', None),
        "primarySpecies" : [getattr(PICCTS_input, "primarySpeciesAq", []),
                           getattr(PICCTS_input, "primarySpeciesPha", []),
                           getattr(PICCTS_input, "primarySpeciesSurf", []),
                           getattr(PICCTS_input, "speciationEch", []),
                           getattr(PICCTS_input, "primarySpeciesPhantom", ['O','H','OH','H2O','pH','ph','pe']),
                           getattr(PICCTS_input, "secSpeciesWaterUndecomposed", []),
                           (getattr(PICCTS_input, "primarySpeciesAq", []) + getattr(PICCTS_input, "primarySpeciesPha", []) +
                            getattr(PICCTS_input, "primarySpeciesSurf", [])+getattr(PICCTS_input, "speciationEch", []))
                           ],
        "SI" : getattr(PICCTS_input, "SI", [0 for pha in getattr(PICCTS_input, "primarySpeciesPha", [])]),
        "mineralReversibility" : getattr(PICCTS_input, "mineralReversibility", ["" for pha in getattr(PICCTS_input, "primarySpeciesPha", [])]),
        "cutoffs" : [getattr(PICCTS_input, 'cutoffPhreeqCAq', 1e-20),
                    getattr(PICCTS_input, 'cutoffPhreeqCPha', -1e-99),
                    getattr(PICCTS_input, 'cutoffPhreeqCSurf', 1e-20),
                    getattr(PICCTS_input, 'cutoffPhreeqCEch', 1e-20)],
        "database": ["",getattr(PICCTS_input, 'masterSpecies', None),
                     getattr(PICCTS_input, 'solutionSpecies', None),
                     phases,
                     getattr(PICCTS_input, 'masterExchange', None),
                     getattr(PICCTS_input, 'exchangeSpecies', None),
                     getattr(PICCTS_input, 'masterSurface', None),
                     getattr(PICCTS_input, 'surfaceSpecies', None),
                     getattr(PICCTS_input, 'kineticRates', None),
                     getattr(PICCTS_input, 'kinetics', None),],
        "fixpH": getattr(PICCTS_input, "fixpH", None),
        "userVarBool": {'pH':getattr(PICCTS_input, 'pH', False),
                        'pe':getattr(PICCTS_input, 'pe', False),
                        'reaction':getattr(PICCTS_input, 'reaction', False),
                        'temperature':getattr(PICCTS_input, 'temperature', False),
                        'alkalinity':getattr(PICCTS_input, 'alkalinity', False),
                        'ionicStrength':getattr(PICCTS_input, 'ionicStrength', False),
                        'water':getattr(PICCTS_input, 'water', False),
                        'charge':getattr(PICCTS_input, 'charge', False),
                        'pourcentError':getattr(PICCTS_input, 'pourcentError', False)},
        "userVarList": {'totals':getattr(PICCTS_input, 'totals', False),
                       'activities':getattr(PICCTS_input, 'activities', False),
                       'saturationIndices':getattr(PICCTS_input, 'saturationIndices', False),
                       'gases':getattr(PICCTS_input, 'gases', False),
                       'kineticReactant':getattr(PICCTS_input, 'kineticReactant', False),
                       'solidSolution':getattr(PICCTS_input, 'solidSolution', False),
                       'isotopes':getattr(PICCTS_input, 'isotopes', False),
                       'calculateValues':getattr(PICCTS_input, 'calculateValues', False)},
        "maillesChargeGeom" : maillesChargeGeom,
        "solMod" : getattr(PICCTS_input, 'solMod', False),
        "speciesCharge" : getattr(PICCTS_input, "speciesCharge", ['pH']),
        "speciesChargeGeometry" : getattr(PICCTS_input, 'speciesChargeGeometry', None),
        "speciationCharge" : getattr(PICCTS_input, "speciationCharge", False),
        "pHdefault" : getattr(PICCTS_input, "pHdefault", 7),
        "tempDefault" : getattr(PICCTS_input, "tempDefault", 25),
        "distAnodCathTotale" : getattr(PICCTS_input, "distAnodCathTotale", None),
        "catholyte" : getattr(PICCTS_input, "catholyte", None),
        "anolyte" : getattr(PICCTS_input, "anolyte", None),
        "solModCharge" : getattr(PICCTS_input, "solModCharge", 0),
        "renouvellement" : getattr(PICCTS_input, 'renouvellement', None),
        "acidicTrspt": getattr(PICCTS_input, 'acidicTrspt', False),
        "PhreeqCClockTime": 0,
        "PhreeqCPrcsTime": 0,
        "water" : getattr(PICCTS_input, 'water', False),
        "kinetics" : getattr(PICCTS_input, 'kinetics', False),
        'supplementarySolution' : getattr(PICCTS_input, 'supplementarySolution', None),
        })
        
        primToSecSpecies = {}
        for _, row in centralDict['commMtrx'][centralDict['systemSpeciation']].iterrows():
            for comp, conc in row.items():
                if isinstance(conc, str):
                    print("\nTypeError: '<' not supported between instances of 'str' and 'int'")
                    print(conc,comp)
                    sys.exit()
                composition, _ = decomposingIntoPrimSpecies(comp, centralDict['primarySpecies'][-1])
                primToSecSpecies.update({comp : composition})

        centralDict.update({"primToSecSpecies" : primToSecSpecies,})
        
    elif centralDict["couplingInfo"][1] == 'xGEMS':
        import gems
        speciationLauncher = {
            2: gems.spct
        }
        centralDict.update({
        "chemPath" : getattr(PICCTS_input, 'chemPath', 'dat.lst'),
        "independentComponents" : getattr(PICCTS_input, "independentComponents", None),
        "xGEMSClockTime" : 0,
        "constantSpecies" : getattr(PICCTS_input, 'constantSpecies', {}),
        "xGEMSPrcsTime" : 0,})
        
        primToSecSpecies = {}
        for _, row in centralDict['commMtrx'][centralDict['systemSpeciation']].iterrows():
            for comp, conc in row.items():
                if isinstance(conc, str):
                    print("\nTypeError: '<' not supported between instances of 'str' and 'int'")
                    print(conc,comp)
                    sys.exit()
                composition, _ = decomposingIntoPrimSpecies(comp, centralDict['independentComponents'])
                primToSecSpecies.update({comp : composition})
        
        userprimToSecSpecies = getattr(PICCTS_input, 'nonTrivialDC', None)
        if userprimToSecSpecies:
            for ky in primToSecSpecies:
                for subky in userprimToSecSpecies:
                    if subky == ky:
                          primToSecSpecies[ky] = userprimToSecSpecies[ky]  
        centralDict.update({"nonTrivialDC" : primToSecSpecies})
        
    elif centralDict["couplingInfo"][1] == 'ORCHESTRA':
        import orchestra
        speciationLauncher = {
            3: orchestra.spct
        }
        
        speciesAttributes = {}
        dico = getattr(PICCTS_input, 'speciesAttributes', { })
        for ky in dico.keys():
            if isinstance(dico[ky], dict) : 
                for subkey in dico[ky]:
                    speciesAttributes[subkey] = f"{dico[ky][subkey]}.{ky}"
                    
            else:
                for spc in centralDict['systemSpeciation']:
                    if spc in dico[ky]:
                        speciesAttributes[spc] = f"{spc}.{ky}"
                        continue
        
        centralDict.update({
        "chemPath" : getattr(PICCTS_input, 'chemPath', 'chemistry1.inp'),
        "primarySpecies" : [getattr(PICCTS_input, "primarySpeciesAq", []),
                           getattr(PICCTS_input, "primarySpeciesPha", []),
                           getattr(PICCTS_input, "primarySpeciesSurf", []),
                           getattr(PICCTS_input, "speciationEch", []),
                           getattr(PICCTS_input, "primarySpeciesPhantom", ['OH','H']),
                           (getattr(PICCTS_input, "primarySpeciesAq", []) + getattr(PICCTS_input, "primarySpeciesPha", []) +
                            getattr(PICCTS_input, "primarySpeciesSurf", [])+getattr(PICCTS_input, "speciationEch", []))
                           ],

        "primarySpeciesEchSorbed" :  getattr(PICCTS_input, "primarySpeciesEchSorbed", {}),
        "ORCHESTRAClockTime": 0,
        "ORCHESTRAPrcsTime": 0,
        })
        
        centralDict.update({
        "inputVariableOrchestra" : [f"{spc}.tot" for spc in centralDict['primarySpecies'][-1]],
        "outputVariableOrchestra" : [f"{speciesAttributes[spc]}" for spc in centralDict['systemSpeciation']],
            })

        primToSecSpecies = {}
        for _, row in centralDict['commMtrx'][centralDict['systemSpeciation']].iterrows():
            for comp, conc in row.items():
                if isinstance(conc, str):
                    print("\nTypeError: '<' not supported between instances of 'str' and 'int'")
                    print(conc,comp)
                    sys.exit()
                composition, _ = decomposingIntoPrimSpecies(comp, centralDict['primarySpecies'][-1])
                primToSecSpecies.update({comp : composition})

        centralDict.update({"primToSecSpecies" : primToSecSpecies,})

    else: 
        print(f'No module is associated with chemModule = {getattr(PICCTS_input, "chemModule", 1)}')
        sys.exit()


    if centralDict["couplingInfo"][2] == 'COMSOL':
        import comsol
        transportLauncher = {
            1: comsol.trspt, 
        }
        centralDict.update({
        "comsolTags" : [getattr(PICCTS_input, 'comsolComp', 'comp1'),
                       getattr(PICCTS_input, 'comsolIntFonction', 'int1'),
                       getattr(PICCTS_input, "comsolStudy", 'std1'), 
                    ],
        "comsolCoupling": getattr(PICCTS_input, 'comsolCoupling', 'data1'),
        "comsolData": getattr(PICCTS_input, 'comsolData', ['data1']),
        "comsolVTU" : getattr(PICCTS_input, 'comsolVTU', True),
        "comsolCore" : getattr(PICCTS_input, 'comsolCore', None),
        "inputPath" : inputPath,
        "COMSOLClockTime" : 0,
        "COMSOLopeningClosing" : 0,
        })
    else:
        print(f'No module is associated with trsptModule = {getattr(PICCTS_input, "trsptModule", 1)}')
        sys.exit()

    print(f"""
####   #  ####  ####  #####  ####     {centralDict['couplingInfo'][0]} splitting :
#  #   #  #     #       #    #        {centralDict['couplingInfo'][1]}--->
####   #  #     #       #    ####            <---{centralDict['couplingInfo'][2]}
#      #  #     #       #       #     {maxTime-timeStepReprise}{centralDict['timeUnit']} in {len(dtPICCTS)-stepReprise} steps        
#      #  ####  ####    #    ####     
       """, flush=True)

    
    if not stepReprise:
        with open("warning.log", "w") as warningLog:
            warningLog.write(f"PICCTS : {centralDict['couplingInfo'][0]} {centralDict['couplingInfo'][1]}-{centralDict['couplingInfo'][2]}.\n")
            if getattr(PICCTS_input, "description", None): warningLog.write(f"PICCTS : Description :\n{PICCTS_input.description}\n")
            warningLog.write(f"PICCTS : Time steps ({centralDict['timeUnit']}) :\n")    
            for dtt in dtPICCTS:
                if dtt == dtPICCTS[-1]: warningLog.write(f"{dtt}.\n")
                else : warningLog.write(f"{dtt}, ")
    else:
        with open("warning.log", "a") as warningLog:
            warningLog.write(f"PICCTS, time = {dtPICCTS[stepReprise]}{centralDict['timeUnit']} : run reprise ...\n")


    for l,t in enumerate(dtPICCTS[start:], start = startTimeStep):
        # print(t)
        if any('time' in x for x in (centralDict['partialCoupledParameters'],centralDict['coupledParametersNames'])):
            centralDict['commMtrx']['time'] = t
        
        pd.set_option('display.max_columns', None)
        pd.set_option('display.max_rows', 10)
        
        if l==startTimeStep and not stepReprise: dt = t
        else: dt = t-dtPICCTS[l-1]

        centralDict.update({
            "dtStep":dt,
            "tStep":t,
            "lStep":l,
            })
        
        print(f"#######  step n°{l+1}/{len(dtPICCTS)}, dt = {dt}{centralDict['timeUnit']}  #######")

        if centralDict["firstStepEquilibrium"] and l==startTimeStep:
            centralDict["dtStep"]=0
            centralDict.update(speciationLauncher.get(centralDict['chemModule'],0)(centralDict))
            # print(centralDict['commMtrx'])
            # sys.exit()
            centralDict["firstStepEquilibrium"]=False
            centralDict["dtStep"]=dt

        if centralDict['couplingInfo'][0]=='SNIA':
            centralDict.update(transportLauncher.get(centralDict['trsptModule'],0)(centralDict))
            # print(centralDict['commMtrx'])
            centralDict.update(speciationLauncher.get(centralDict['chemModule'],0)(centralDict))    
            # print(centralDict['commMtrx'])
        elif centralDict['couplingInfo'][0]=='Strang':
            centralDict['dtStep'] = dt/2
            centralDict.update(transportLauncher.get(centralDict['trsptModule'],0)(centralDict))
            centralDict['dtStep'] = dt
            centralDict.update(speciationLauncher.get(centralDict['chemModule'],0)(centralDict))
            centralDict['dtStep'] = dt/2
            centralDict.update(transportLauncher.get(centralDict['trsptModule'],0)(centralDict))

        elif centralDict['couplingInfo'][0]=='Alternative':
            if l%2==0:
                centralDict.update(transportLauncher.get(centralDict['trsptModule'],0)(centralDict))
                centralDict.update(speciationLauncher.get(centralDict['chemModule'],0)(centralDict))
            else:
                centralDict.update(speciationLauncher.get(centralDict['chemModule'],0)(centralDict))
                centralDict.update(transportLauncher.get(centralDict['trsptModule'],0)(centralDict))

        elif centralDict['couplingInfo'][0]=='Additive':
            commMtrxIC = centralDict['commMtrx'].copy()
            centralDict.update(transportLauncher.get(centralDict['trsptModule'],0)(centralDict))
            commMtrxTrspt = centralDict['commMtrx'].copy()
            centralDict['commMtrx'] = commMtrxIC.copy()
            centralDict.update(speciationLauncher.get(centralDict['chemModule'],0)(centralDict))

            commMtrx_Final = centralDict['commMtrx'] + commMtrxTrspt - commMtrxIC 
            centralDict['commMtrx'] = commMtrx_Final.copy()
            
        elif centralDict['couplingInfo'][0]=='Symmetrical': # Output solely for last symmetrical coupling ..
            commMtrxIC = centralDict['commMtrx'].copy()
            centralDict.update(transportLauncher.get(centralDict['trsptModule'],0)(centralDict))
            centralDict.update(speciationLauncher.get(centralDict['chemModule'],0)(centralDict))
            commMtrxSym1 = centralDict['commMtrx'].copy()
            print('\\')
            centralDict['commMtrx'] = commMtrxIC.copy() 
            centralDict.update(speciationLauncher.get(centralDict['chemModule'],0)(centralDict))
            centralDict.update(transportLauncher.get(centralDict['trsptModule'],0)(centralDict))
            
            commMtrx_Final = (commMtrxSym1 + centralDict['commMtrx'])/2
            centralDict['commMtrx'] = commMtrx_Final.copy()
        
        stepTime = time.time() - PICCTSstart
        
        centralDict['commMtrx'].to_csv(os.path.join(paths['outputHistory'], f"CommMtrx_{l+1}.txt"), index=False, header=True, sep='\t')

        if t!=max(dtPICCTS):    
            print(f"Remaining calculation time ~ {writeTime((max(dtPICCTS) * stepTime - timeStepReprise)/ (t-timeStepReprise) - stepTime,2)}")

    print()
    with open("warning.log", "a") as warningLog:
        chem = centralDict["couplingInfo"][1]
        trspt = centralDict["couplingInfo"][2]
        warningLog.write(f"PICCTS : Total calculation time : {writeTime(time.time() - PICCTSstart - centralDict['waitingTime'])}\n")
        warningLog.write(f"PICCTS : {writeTime(centralDict[f'{chem}PrcsTime'])} of processor time spend by {chem}\n")
        warningLog.write(f"PICCTS : {writeTime(centralDict[f'{chem}ClockTime'])} of wall clock time spend by {chem}\n")
        warningLog.write(f"PICCTS : {writeTime(centralDict[f'{trspt}ClockTime'])} spend by {trspt}\n")
        if centralDict["trsptModule"] == 1:
            warningLog.write(f"PICCTS : {writeTime(centralDict[f'{trspt}openingClosing'])} of COMSOL file processing, spend by both {trspt} and PICCTS\n")
            warningLog.write(f"PICCTS : {writeTime(time.time() - PICCTSstart -centralDict[f'{chem}ClockTime'] - centralDict[f'{trspt}ClockTime'] - centralDict['waitingTime'] - centralDict[f'{trspt}openingClosing'])} of interfacing spend by PICCTS\n") # wall clock time
        else:
            warningLog.write(f"PICCTS : {writeTime(time.time() - PICCTSstart -centralDict[f'{chem}ClockTime'] - centralDict[f'{trspt}ClockTime'] - centralDict['waitingTime'])} of interfacing spend by PICCTS\n") # wall clock time

        if centralDict['waitingTime'] : warningLog.write(f"PICCTS : {writeTime(centralDict['waitingTime'])} spend by {trspt}\n")

    if centralDict['warningRun']:
        with open("warning.log", "a") as warningLog: warningLog.write(f"PICCTS : Total warnings :\n{centralDict['warningRun']}")
        print(f"Oh no, {centralDict['warningRun']} warnings occured :( ! See warning.log file.")

    if getattr(PICCTS_input, "description", None): print(f"Description = {PICCTS_input.description}")