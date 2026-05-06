import numpy as np
import time 
import pandas as pd
import sys 
import importlib.util
import os
import concurrent.futures

import PyORCHESTRA

from pathlib import Path
current_dir = Path(__file__).parent
if str(current_dir) not in sys.path:
    sys.path.insert(0, str(current_dir))

with open("store.txt", 'r', encoding='utf-8') as fichier:
    inputPath =  fichier.readline().strip() 
    nameInput =  fichier.readline().strip()

PICCTSstart = time.time()

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

def speciationOrchestra(orchDict,commMtrxPart):
    """
    Site quantity is natively defined in the .inp file. This file could change depending on the node number.
    side note:
        Sorbed attribute is .solid
        Speciation is Na+.con, Na.diss en solution, Na.tot total. Phase attributes are defined by speciesAttributes keyword. They may refer to the .inp file.
    """
    
    commMtrx_species = commMtrxPart[orchDict['systemSpeciation']].copy()
    
    
    dico = {spc: [0]*len(orchDict['commMtrx']) for spc in orchDict['primarySpecies'][-1]}  
    line = 0
    for _, row in commMtrx_species.iterrows():
        for comp, conc in row.items():
            for prim in orchDict['primToSecSpecies'][comp] :
                dico[prim][line] += orchDict['primToSecSpecies'][comp][prim] * conc
        line += 1

    commMtrx_primSpecies = pd.DataFrame(dico)

    solver = PyORCHESTRA.ORCHESTRA()
    solver.initialise(orchDict['chemPath'], 1, orchDict['inputVariableOrchestra'], orchDict['outputVariableOrchestra'])
        
    outputOrchestra = pd.DataFrame(columns = orchDict['outputVariableOrchestra'])
    rows = commMtrx_primSpecies.to_numpy()
    
    results = []
    calcTime = 0
    
    for row in rows:
        ref = time.perf_counter()
        results.append(solver.set_and_calculate_single(row)[0])
        calcTime += time.perf_counter() - ref
    
    outputOrchestra = pd.DataFrame(results, columns=orchDict['outputVariableOrchestra'])
    outputOrchestra.columns = orchDict['systemSpeciation']
    outputOrchestra = pd.concat([commMtrxPart[orchDict["coord"]],outputOrchestra], axis=1)

    return outputOrchestra, commMtrx_primSpecies, calcTime



def spct(centralDict):
    print("ORCHESTRA", end=" ", flush=True)
    startOrchestra = time.time()
    chunk_size = int(np.ceil(len(centralDict["commMtrx"]) / centralDict["PIDnbr"]))

    commMtrxSplit = [centralDict["commMtrx"].iloc[i:i + chunk_size] for i in range(0, len(centralDict["commMtrx"]), chunk_size)]
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=centralDict["PIDnbr"]) as executor:
        futures = []
        for chunk in commMtrxSplit:
            futures.append(
                executor.submit(
                    speciationOrchestra,
                    centralDict,
                    chunk
                )
            )
            
    results = []
    
    results = [f.result() for f in futures]

    df1, df2, intgr2 = zip(*results)

    commMtrxSpct = pd.concat(df1, ignore_index=True)
    commMtrx_primSpecies = pd.concat(df2, ignore_index=True)
    prcsTime = sum(intgr2)
    
    centralDict.update({
        "commMtrx": commMtrxSpct,
        "ORCHESTRAClockTime": centralDict['ORCHESTRAClockTime'] + time.time() - startOrchestra,
        "ORCHESTRAPrcsTime": centralDict['ORCHESTRAPrcsTime'] + prcsTime,
        })
    
    
    if centralDict["firstStepEquilibrium"]==True:
        commMtrx_primSpecies.to_csv(os.path.join(centralDict['paths']['primSpecies'], f"PrimarySpecies_{centralDict['lStep']}.txt"), index=False, header=True, sep='\t')
        commMtrxSpct.to_csv(os.path.join(centralDict['paths']['outputSpeciation'], f"ORCHESTRA_{centralDict['lStep']}.txt"), index=False, header=True, sep='\t')
    else:
        commMtrx_primSpecies.to_csv(os.path.join(centralDict['paths']['primSpecies'], f"PrimarySpecies_{centralDict['lStep']+1}.txt"), index=False, header=True, sep='\t')
        commMtrxSpct.to_csv(os.path.join(centralDict['paths']['outputSpeciation'], f"ORCHESTRA_{centralDict['lStep']+1}.txt"), index=False, header=True, sep='\t')
   
    print(f"({writeTime((time.time() - startOrchestra))})") 

    return centralDict