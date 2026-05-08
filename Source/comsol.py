import time 
import pandas as pd
import sys 
import importlib.util
import os

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

def transportComsol(comsolDict,queue):

    import mph
    
    scriptWarning = ""
    warning = 0
    abort = False
    calcTime = 0
    COMSOLopeningClosing = 0
    inputPath = os.path.join(comsolDict['inputPath'], 'inputTransport.txt')
    
    comsolDict['commMtrx'].to_csv(inputPath, index=False, header=False, sep='\t') 
    # sys.exit()
    client = mph.start()
    if comsolDict['comsolCore'] : client = mph.Client(cores=comsolDict['comsolCore'])
    
    licenceTaken = False
    licence = False     
    while not licence:
        try: 
            ref = time.perf_counter()
            model = client.load(f"{comsolDict['paths']['trsptPath']}")
            javamodel = model.java
            javamodel.component(f"{comsolDict['comsolTags'][0]}").func(f"{comsolDict['comsolTags'][1]}").set("filename", f"{inputPath}")
            COMSOLopeningClosing += time.perf_counter() - ref
            licence = True
            if licenceTaken:
                endWaiting = time.time()
                print("Licence released :)", end=" ", flush=True)

        except Exception as r :
            if not licenceTaken:
                startWaiting = time.time()
                print(r)
                print("Waiting for a licence ? ...", flush=True)
                licenceTaken = True
                scriptWarning+= f"COMSOL, time = {comsolDict['tStep']}{comsolDict['timeUnit']}, time-step n° {comsolDict['lStep']+1} : {r}\n"
                time.sleep(3)
            if licenceTaken: time.sleep(3)

    if licenceTaken:
        warning += 1
        waitingLicence = endWaiting - startWaiting
        print(f"({writeTime(waitingLicence)} of waiting)", flush = True)
        scriptWarning+=f"COMSOL, t={comsolDict['tStep']}{comsolDict['timeUnit']}, time-step n° {comsolDict['lStep']+1} : {writeTime(waitingLicence)} of waiting \n"
        waitingLicence = endWaiting - startWaiting
    else: waitingLicence=0 

    # some Java to dynamically couple COMSOL ..
    # COMSOL output should fit the expected communication matrix (cf user manuel)
    try:
        ref = time.perf_counter()
        javamodel.component(f"{comsolDict['comsolTags'][0]}").func(f"{comsolDict['comsolTags'][1]}").refresh();
        COMSOLopeningClosing += time.perf_counter() - ref
        if comsolDict['timeUnit'] =='y': timeUnit = 'a'
        else: timeUnit = comsolDict['timeUnit']
        ref = time.perf_counter()
        javamodel.study(f"{comsolDict['comsolTags'][2]}").feature("time").set("tunit", f"{timeUnit}");
        javamodel.study(f"{comsolDict['comsolTags'][2]}").feature("time").set("tlist", f"range(0,{comsolDict['dtStep']},{comsolDict['dtStep']})");
        COMSOLopeningClosing += time.perf_counter() - ref
    except Exception as r:
        print(r)
        warning += 1
        scriptWarning+= f"COMSOL, time = {comsolDict['tStep']}{comsolDict['timeUnit']}, time-step n° {comsolDict['lStep']+1} : Fatal error occured : \n{r}\n"
        abort = True
    
    if not abort:
        try:
            ref = time.perf_counter()
            model.solve()
            calcTime += time.perf_counter() - ref
        except Exception as r:
            print(r)
            with open("warning.log", "a") as warningLog: warningLog.write(f"COMSOL, t={comsolDict['tStep']}{comsolDict['timeUnit']}, time-step n° {comsolDict['lStep']} : {r}\n")
            abort = True
    if not abort:
        outputPath = os.path.join(comsolDict['inputPath'], 'outputTransport.txt')
        ref = time.perf_counter()
        javamodel.result().export(f"{comsolDict['comsolCoupling']}").setIndex("looplevelinput", "last", 0);
        javamodel.result().export(f"{comsolDict['comsolCoupling']}").set("exporttype", "text");
        javamodel.result().export(f"{comsolDict['comsolCoupling']}").set("filename", f"{outputPath}");
        javamodel.result().export(f"{comsolDict['comsolCoupling']}").run();
        COMSOLopeningClosing += time.perf_counter() - ref
        
        for i,tag in enumerate(comsolDict['comsolData']): # user defined variables ..
            ref = time.perf_counter()
            path = os.path.join(comsolDict['paths'][f'outputTransport{tag}'], f'COMSOL_{comsolDict["lStep"]+1}.txt')
            javamodel.result().export(f"{tag}").setIndex("looplevelinput", "last", 0);
            javamodel.result().export(f"{tag}").set("exporttype", "text");
            javamodel.result().export(f"{tag}").set("filename", f"{path}");
            javamodel.result().export(f"{tag}").run();
            COMSOLopeningClosing += time.perf_counter() - ref
            if comsolDict['comsolVTU']:
                path = os.path.join(comsolDict['paths'][f'VTU{tag}'],f'COMSOL_{comsolDict["lStep"]+1}.vtu')
                ref = time.perf_counter()
                javamodel.result().export(f"{tag}").setIndex("looplevelinput", "last", 0);
                javamodel.result().export(f"{tag}").set("exporttype", "vtu");
                javamodel.result().export(f"{tag}").set("filename", f"{path}");
                javamodel.result().export(f"{tag}").run();
                COMSOLopeningClosing += time.perf_counter() - ref
                if comsolDict['lStep'] == 0:
                    path = os.path.join(comsolDict['paths'][f'VTU{tag}'],'COMSOL_0.vtu')
                    ref = time.perf_counter()
                    javamodel.result().export(f"{tag}").setIndex("looplevelinput", "first", 0);
                    javamodel.result().export(f"{tag}").set("exporttype", "vtu");
                    javamodel.result().export(f"{tag}").set("filename", f"{path}");
                    javamodel.result().export(f"{tag}").run();
                    COMSOLopeningClosing += time.perf_counter() - ref

        ref = time.perf_counter()
        client.clear()
        COMSOLopeningClosing += time.perf_counter() - ref
        
        commMtrx_Comsol = pd.read_csv(outputPath,sep=r"\s+",comment="%",header=None, names=list(comsolDict['commMtrx'].columns))
        results = [commMtrx_Comsol, warning, scriptWarning,waitingLicence, abort,calcTime]
        
        if queue : queue.put(results)
        return commMtrx_Comsol, warning, scriptWarning,waitingLicence, abort,calcTime, COMSOLopeningClosing

    else:
        return pd.DataFrame(), warning, scriptWarning, waitingLicence, abort, 0, COMSOLopeningClosing

def trspt(centralDict):
    print("COMSOL", end=" ", flush=True)
    startComsol = time.time()
    
    comm,  warning, scriptWarning,waitingLicence, abort, calcTime, COMSOLopeningClosing = transportComsol(
        centralDict, None)
    
    if warning:
        with open("warning.log", "a") as warningLog:
            warningLog.write(f"COMSOL, time = {centralDict['tStep']}{centralDict['timeUnit']}, time step n°{centralDict['lStep']+1} : the {warning} following warnings occured ...\n {scriptWarning}")
    if abort:
        with open("warning.log", "a") as warningLog:
            warningLog.write(f"COMSOL, time = {centralDict['tStep']}{centralDict['timeUnit']}, time step n°{centralDict['lStep']+1} : Fatal COMSOL error. Aborting run.")
        print('Fatal COMSOL error. Aborting run.')
        sys.exit()

    centralDict.update({
            "commMtrx": comm,
            "COMSOLClockTime": centralDict['COMSOLClockTime'] + calcTime,
            "waitingTime" : centralDict['waitingTime'] + waitingLicence,
            "COMSOLopeningClosing" : centralDict['COMSOLopeningClosing'] + COMSOLopeningClosing,
            })

    print(f"({writeTime((time.time()-startComsol))})") 

    return centralDict
