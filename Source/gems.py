import xgems
import pandas as pd
import sys
import numpy as np
import time 
import importlib.util
import os
import concurrent.futures

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

def writeTime(tps, arr=2):
    if tps >= 3600 * 24:
        return f"{tps / (3600 * 24):.{arr}f} d"
    elif tps >= 3600:
        return f"{tps / 3600:.{arr}f} h"
    elif tps >= 60:
        return f"{tps / 60:.{arr}f} min"
    else:
        return f"{tps:.{arr}f} sec"

gemsStatus= {
0: "No GEM re-calculation needed",
1: "Need GEM calculation with LPP (automatic) initial approximation (AIA)",
2: "OK after GEM calculation with LPP AIA",
3: "Bad (not fully trustful) result after GEM calculation with LPP AIA",
4: "Failure (no result) in GEM calculation with LPP AIA",
5: "Need GEM calculation with no-LPP (smart) IA, SIA using the previous speciation",
6: "OK after GEM calculation with SIA",
7: "Bad (not fully trustful) result after GEM calculation with SIA",
8: "Failure (no result) in GEM calculation with SIA",
9: "Terminal error in GEMS3K (e.g., memory corruption). Restart required.",
    }

def speciation_xGEMS(gemsDict,commMtrxPart):

    engine = xgems.ChemicalEngine(gemsDict['chemPath'])
    commMtrxCoord = commMtrxPart[gemsDict['coord']].copy()
    commMtrx_species = commMtrxPart.drop(columns=gemsDict['coord']) # To be sure to not change nodes coordinates
    outputGems = pd.DataFrame(columns=gemsDict['systemSpeciation'])

    dico = {}
    dico = {spc: [0]*len(commMtrx_species) for spc in gemsDict['independentComponents']}


    ligne = 0
    for _, row in commMtrx_species.iterrows():
        for comp, conc in row.items():
            for prim in gemsDict['nonTrivialDC'][comp] :
                dico[prim][ligne] += gemsDict['nonTrivialDC'][comp][prim] * conc
        ligne += 1

    commMtrx_primSpecies = pd.DataFrame(dico)
    commMtrx_primSpecies.index = commMtrx_species.index
    
    commMtrx_primSpecies = commMtrx_primSpecies.astype(float)
    commMtrx_primSpecies= commMtrx_primSpecies.clip(lower=1e-16)

    gemsStatusList = []
    gemsIterations = []
        
    calcTime = 0
    for index in commMtrx_primSpecies.index:
        ref = time.perf_counter()
        status = engine.equilibrate(
            298,1e5, # to continue ..
            [
                gemsDict['constantSpecies'][spc] if spc in gemsDict['constantSpecies']
                else commMtrx_primSpecies.at[index, spc]
                for spc in gemsDict['independentComponents']
            ])
        calcTime += time.perf_counter() - ref
            
        outputGems.loc[len(outputGems)] = engine.speciesAmounts()
        gemsStatusList += [status]
        gemsIterations += [engine.numIterations()]
        
    outputGems.index = commMtrx_species.index
    outputGems = pd.concat([commMtrxCoord,outputGems], axis=1)
    
    return outputGems,commMtrx_primSpecies,gemsStatusList,gemsIterations, calcTime
	
def spct(centralDict): 
    print("xGEMS", end=" ", flush=True)
    startGems = time.time()
    chunk_size = int(np.ceil(len(centralDict['commMtrx']) / centralDict['PIDnbr']))

    commMtrxSplit = [centralDict['commMtrx'].iloc[i:i + chunk_size] for i in range(0, len(centralDict['commMtrx']), chunk_size)]
    with concurrent.futures.ProcessPoolExecutor(max_workers=centralDict['PIDnbr']) as executor:
        futures = []
        for chunk in commMtrxSplit:
            futures.append(
                executor.submit(
                    speciation_xGEMS,
                    centralDict,
                    chunk,
                )
            )
    
   
    results = []   
    
    results = [f.result() for f in futures]

    df1,df2, status, iteration, intgr2= zip(*results)

    commMtrxPart = pd.concat(df1, ignore_index=True)
    commMtrx_primSpecies = pd.concat(df2, ignore_index=True)
    prcsTime = sum(intgr2)
    
    gemsStatusList = status[0]
    gemsIterations = iteration[0]


    if centralDict["firstStepEquilibrium"]==True:
        commMtrx_primSpecies.to_csv(os.path.join(centralDict['paths']['primSpecies'], f"PrimarySpecies_{centralDict['lStep']}.txt"), index=False, header=True, sep='\t')
        commMtrxPart.to_csv(os.path.join(centralDict['paths']['outputSpeciation'], f"PhreeqC_{centralDict['lStep']}.txt"), index=False, header=True, sep='\t')
    else:
        commMtrx_primSpecies.to_csv(os.path.join(centralDict['paths']['primSpecies'], f"PrimarySpecies_{centralDict['lStep']+1}.txt"), index=False, header=True, sep='\t')
        commMtrxPart.to_csv(os.path.join(centralDict['paths']['outputSpeciation'], f"PhreeqC_{centralDict['lStep']+1}.txt"), index=False, header=True, sep='\t')
       
    with open("warning.log", "a") as warningLog:
        for i,status in enumerate(gemsStatusList):
            if status!=2:
                warningLog.write(f"GEMS, node n°{i}, time step n°{centralDict['lStep']+1}, t={centralDict['tStep']}{centralDict['timeUnit']}, PID={os.getpid()} : with {gemsIterations[i]} :\n")
                warningLog.write(f"GEMS, node n°{i}, time step n°{centralDict['lStep']+1}, t={centralDict['tStep']}{centralDict['timeUnit']}, PID={os.getpid()} : {gemsStatus[status]}\n")

    
    if centralDict["firstStepEquilibrium"]==True:
        commMtrx_primSpecies.to_csv(os.path.join(centralDict['paths']['primSpecies'], f"PrimarySpecies_{centralDict['lStep']}.txt"), index=False, header=True, sep='\t')
    else:
        commMtrx_primSpecies.to_csv(os.path.join(centralDict['paths']['primSpecies'], f"PrimarySpecies_{centralDict['lStep']+1}.txt"), index=False, header=True, sep='\t')

    print(f"({writeTime((time.time() - startGems))})") 
    
    centralDict.update({
        "commMtrx": commMtrxPart,
        "xGEMSClockTime": centralDict['xGEMSClockTime'] + time.time() - startGems,
        "xGEMSPrcsTime": centralDict['xGEMSPrcsTime'] + prcsTime,
        })
 

    return centralDict