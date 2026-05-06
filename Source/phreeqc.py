import numpy as np
import time 
import pandas as pd
import sys 
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


def speciesToNode(v, ranges, labels):
    """
    Associate a species to a node number
    """
    ranges = np.asarray(ranges, dtype=float)  
    
    results = []
    
    mask = (v >= ranges[:, 0]) & (v <= ranges[:, 1])
    
    if np.any(mask):
        idx = np.where(mask)[0][0]
        matched_label =  labels[idx]
    else:
        matched_label = None
    
    results.append(matched_label)
    
    return results


def speciationPhreeqC(pcDict, commMtrxPart,rnvlmt): 
    
    import phreeqpy.iphreeqc.phreeqc_dll as phreeqc_mod
    phreeqc = phreeqc_mod.IPhreeqc()
    phreeqc.load_database(pcDict['paths']['chemPath'])
    
    systemSpeciation = [cell for cell in list(commMtrxPart.columns) if cell not in ['x','y','z']] # répétitif ...
    
    def nodeSpeciation(commMtrxPart,commMtrx_primSpecies,spcChargeDefaut,spcChargeGeom,speciationCharge,solMod): # complexationSurface, echangeIon, phases,primarySpeciesSurf,primarySpeciesAq, primarySpeciesAq, primarySpeciesPha en variables 'globales'

        phListUser = ['ph','pH','H']
        scriptMain =f"node n°{index}\n"
        if pcDict['water']:
            for key in pcDict['water']:
                if index <= key:
                    water = pcDict['water'][key]
                    break
        else: water = 1
        
        scriptAnodeCell2 ="""\nSOLUTION 1
    \t-units mol/kgw
    pH 2 charge
    S 0.005
    Ammo 0.09
    Cl 0.09
    Traceur 0.1\n"""    

        scriptCathodeCell2 ="""\nSOLUTION 1
    \t-units mol/kgw
    pH 7 charge
    Acetate 0.1
    Ammo 0.1\n"""  


        scriptAnodeCell1 ="""\nSOLUTION 1
    \t-units mmol/kgw
    pH 7 charge
    Ca 0.6
    Mg 0.15
    Na 3.1
    K 0.05
    Cl 3.15
    S 0.15
    C 1.2
    Traceur 0.1\n"""    

        scriptCathodeCell1 ="""\nSOLUTION 1
    \t-units mmol/kgw
    pH 7 charge
    Ca 0.6
    Mg 0.15
    Na 3.1
    K 0.05
    Cl 3.15
    S 0.15
    C 1.2\n"""  
        
        bait = f"""
    SOLUTION 1
    -units mmol/kgw
    # Na 1
    # Cl 1
    -water {water}
    END
    
    RUN_CELLS
    -cells 1
    END               
        """
        if pcDict['kinetics'] :
            scriptMain += f"KINETICS\n{pcDict['kinetics']}"
            scriptMain += f'\n-step {pcDict["dtStep"]}\n' #kinetics

        if pcDict['fixpH'] or pcDict['primarySpecies'][1]:
            scriptMain +="\nEQUILIBRIUM_PHASES 1\n"
            if pcDict['fixpH']:
                scriptMain += f"Fix_ph {pcDict['fixpH']}\n"
            if pcDict['primarySpecies'][1] : 
                for j,phase in enumerate(pcDict['primarySpecies'][1]):
                    if commMtrx_primSpecies.loc[index ,phase] > pcDict['cutoffs'][1]:
                        scriptMain +=f"\t{phase} {pcDict['SI'][j]} {commMtrx_primSpecies.loc[index,phase]}  {pcDict['mineralReversibility'][j]}\n"
                # if pcDict['system'] == 'sebCell1' or pcDict['system'] == 'sebCell2' :
                    # if 'Calcite' in systemSpeciation:
                        # scriptMain += f"\tCalcite 0 {commMtrx_primSpecies.loc[index,'Calcite']} dissolve_only\n"
                    # if  index > pcDict['distAnodCathTotale']/(5e-4) and index < ( pcDict['distAnodCathTotale']+ 15e-2)/(5e-4): scriptMain += "\tCO2(g) -3.5 1\n"
                
        if pcDict['primarySpecies'][2]:
            scriptMain +="\nSURFACE 1\n"
            edl = False
            for surface in pcDict['primarySpecies'][2]:
                if commMtrx_primSpecies.loc[index,surface] > pcDict['cutoffs'][2] :
                    scriptMain +=f"\t{surface} {commMtrx_primSpecies.loc[index,surface]} 100 1 \n"
                    edl = True
            if edl and "DebyeLength(m)" in list(commMtrxPart.columns):
                scriptMain +=f"-diffuse_layer {commMtrxPart.loc[index,'DebyeLength(m)']}\n"
        
        if pcDict['primarySpecies'][3]:
            scriptMain +="\nEXCHANGE 1\n"
            for exchange in pcDict['primarySpecies'][3]:
                if commMtrx_primSpecies.loc[index,exchange] > pcDict['cutoffs'][3]:
                        scriptMain +=f"\t{exchange} {commMtrx_primSpecies.loc[index,exchange]}\n"
                        
            if pcDict['acidicTrspt'] and not pcDict['rnvllmt']:
                if AciditeDiff.loc[index,'H+'] >0: scriptMain += f"\tAcide_in_H {AciditeDiff.loc[index,'H+']}\n"
                elif AciditeDiff.loc[index,'H+'] <0: scriptMain += f"\tAcide_out_Similication {-AciditeDiff.loc[index,'H+']}\n"
                if AciditeDiff.loc[index,'OH-'] >0: scriptMain += f"\tBase_in_OH {AciditeDiff.loc[index,'OH-']}\n"
                elif AciditeDiff.loc[index,'OH-'] <0: scriptMain += f"\tBase_out_Similianion {-AciditeDiff.loc[index,'OH-']}\n"
         
        if rnvlmt and index <= pcDict['anolyte']:
            if pcDict['system'] =='sebCell1' : scriptMain += scriptAnodeCell1
            if pcDict['system'] =='sebCell2' : scriptMain += scriptAnodeCell2
            solMod = False
            
        elif rnvlmt and index >= pcDict['catholyte']:
            if pcDict['system'] =='sebCell1' : scriptMain += scriptCathodeCell1
            if pcDict['system'] =='sebCell2' : scriptMain += scriptCathodeCell2
        
        elif solMod:
            # if 'H2O' not in commMtrxPart.columns:
                # waterAmount = 55.509
            # else: waterAmount = commMtrxPart.loc[index,'H2O']
            
            scriptMain += bait + "\nSOLUTION_MODIFY 1\n"
            if 'pH' in systemSpeciation and 'pe' in systemSpeciation:
                scriptMain +=f'''
pH {commMtrxPart.loc[index,'pH']}
pe {commMtrxPart.loc[index,'pe']}\n'''
            
            scriptMain += f"-total_h {water * (commMtrx_primSpecies.loc[index,'H'] + commMtrxPart.loc[index,'H2O']*2) }\n"
            scriptMain += f"-total_o {water * (commMtrx_primSpecies.loc[index,'O'] + commMtrxPart.loc[index,'H2O']) }\n"
            
            scriptMain += "-cb 0\n-totals\n"
            for species in pcDict['primarySpecies'][0]:
                if commMtrx_primSpecies.loc[index,species] > pcDict['cutoffs'][0] and species not in pcDict['primarySpecies'][4] and species != 'Potential':
                    scriptMain += f"\t{species} {water * commMtrx_primSpecies.loc[index,species]}\n"
            scriptMain +="\nRUN_CELLS\n-cells 1\n"

        elif not rnvlmt:
            scriptMain +=f"""\nSOLUTION 1 #node n°{index}
-units mol/kgw
-water {water}
temp {pcDict['tempDefault']}\n"""
            if 'pH' in systemSpeciation and 'pe' in systemSpeciation:
                scriptMain +=f'''
pH {commMtrxPart.loc[index,'pH']}
pe {commMtrxPart.loc[index,'pe']}\n'''

            elif pcDict["acidicTrspt"]:
                scriptMain += f"\tpH {AcidicEcho.loc[index,'pH']}"
            else:
                scriptMain += f"\tpH {pcDict['pHdefault']}" 
            if speciationCharge:
                if spcChargeGeom and spcChargeGeom[0][0] in phListUser : scriptMain += "\tcharge"
                elif spcChargeDefaut and spcChargeDefaut[0] in phListUser : scriptMain += "\tcharge"
            if pcDict["acidicTrspt"] and AciditeDiff.loc[index,'H+'] > 0: scriptMain += f"\tSimilication {AciditeDiff.loc[index,'H+']}\n"
            if pcDict["acidicTrspt"] and AciditeDiff.loc[index,'OH-'] > 0 : scriptMain += f"\tSimilianion {AciditeDiff.loc[index,'OH-']}\n"                            
            scriptMain += "\n"
            for species in pcDict['primarySpecies'][0]:
                if commMtrx_primSpecies.loc[index,species] > pcDict['cutoffs'][0] and species not in pcDict['primarySpecies'][4]:
                    scriptMain += f"\t{species} {commMtrx_primSpecies.loc[index,species]}"
                    if speciationCharge:
                        if spcChargeGeom and spcChargeGeom[0][0] == species : scriptMain += "\tcharge\n"
                        elif spcChargeDefaut and spcChargeDefaut[0] == species : scriptMain += "\tcharge\n"
                        else: scriptMain += "\n"
                    else: scriptMain += "\n"
                    
        if pcDict['supplementarySolution'] : scriptMain += pcDict['supplementarySolution'] + '\n'
                        
        if pcDict['coupledParameters'] :
            scriptMain += "\nUSER_PUNCH\n\t-heading"
            for name in pcDict['coupledParameters'].keys():
                val = pcDict['coupledParameters'][name]
                if len(val) > 1 and isinstance(val[1], list):
                    for spc in val[1]:
                        scriptMain += f"\t{name}{spc}"
                else: scriptMain += f"\t{name}"
            scriptMain += '\n'
            it = 1
            for value in pcDict['coupledParameters'].values():
                if len(value) > 1 and all(isinstance(v, list) for v in value[1:]):
                    for a,mainCpld in enumerate(value[1]):
                        scriptMain += f'\t{it} PUNCH {value[0]}("{mainCpld}"'
                        if len(value) > 2:
                            for subLst in value[2:]:
                                scriptMain += f',"{subLst[a]}"'
                        scriptMain += ")\n"
                        it+=1
                else: 
                    scriptMain += f"\t{it} PUNCH {value[0]}\n"
                    it+=1
        
        scriptMain += """\n\nSELECTED_OUTPUT
    -reset false \n"""
        if 'pH' in systemSpeciation and 'pe' in systemSpeciation:
            scriptMain += """
            -pH True
            -pe True\n"""
        scriptMain +="-molalities"
        for espece in systemSpeciation: 
            if espece not in ['pH','pe','Potential']:    
                if not pcDict['primarySpecies'][1]: scriptMain += f"\t{espece}"
                else:
                    if espece in pcDict['primarySpecies'][1]: pass
                    else: scriptMain += f"\t{espece}"
    
        if pcDict['primarySpecies'][1]:
            scriptMain += "\n\t-equilibrium_phases"
            for phase in pcDict['primarySpecies'][1]: scriptMain +=f"\t{phase}" 
        scriptMain += "\n"
        for val in pcDict['userVarBool'].keys():
            if pcDict['userVarBool'][val]:
                scriptMain += f"-{val} True\n"
        for val in pcDict['userVarList'].keys():
            if pcDict['userVarList'][val]:
                scriptMain += f"-{val}"
                for lst in pcDict['userVarList'][val]:
                    scriptMain += f"\t{lst}"
    
        scriptMain += "\nEND\n"
        # if index == 300:
        #     print(scriptMain)
        #     sys.exit()
        return scriptMain
   
    
    warnings = 0
    scriptWarnings =""
    
    commMtrx_species = commMtrxPart[systemSpeciation].copy()
    # print()
    # print(commMtrxPart)
    # print()
    # print(commMtrx_species)
    # sys.exit()

    dico = {}
    dico = {spc: [0]*len(pcDict['commMtrx']) for spc in pcDict['primarySpecies'][-1]}    


    ligne = 0
    for _, row in commMtrx_species.iterrows():
        for comp, conc in row.items():
            for prim in pcDict['primToSecSpecies'][comp] :
                dico[prim][ligne] += pcDict['primToSecSpecies'][comp][prim] * conc
        ligne += 1
    
    pd.set_option('display.max_rows', None)
    commMtrx_primSpecies = pd.DataFrame(dico)
    # print(commMtrx_primSpecies['Calcite'])
    # sys.exit()
    
    if ['H+'] in systemSpeciation:
        commMtrx_primSpecies['H'] = commMtrx_species['H+']
    if ['OH-'] in systemSpeciation:
        commMtrx_primSpecies['OH'] = commMtrx_species['OH-']

    if pcDict['acidicTrspt']:

        if pcDict['AcidicEcho'].empty:
            pcDict['AcidicEcho'] = pd.DataFrame(columns = ['H+','OH-','pH'], index = commMtrxPart.index)
            pcDict['AcidicEcho']['OH-'] = 0
            pcDict['AcidicEcho']['H+'] = 0
            pcDict['AcidicEcho']['pH'] = -np.log10(commMtrx_species['H+']) # environ.
            AciditeDiff = pcDict['AcidicEcho'].copy()
        else:
            pcDict['AcidicEcho'] = pcDict['AcidicEcho'].loc[commMtrx_species.index]
            AciditeDiff = commMtrx_species[['H+', 'OH-']] - pcDict['AcidicEcho'][['H+', 'OH-']]
        
        AciditeDiff.index = commMtrx_species.index
        
        mask = (AciditeDiff["H+"] > 0) & (AciditeDiff["OH-"] > 0)
        
        min_vals = AciditeDiff.loc[mask, ["H+", "OH-"]].min(axis=1)
        max_vals = AciditeDiff.loc[mask, ["H+", "OH-"]].max(axis=1)
        
        AciditeDiff.loc[mask, "H+"] = np.where(
            AciditeDiff.loc[mask, "H+"] == min_vals, 0, max_vals - min_vals
        )
        
        AciditeDiff.loc[mask, "OH-"] = np.where(
            AciditeDiff.loc[mask, "OH-"] == min_vals, 0, max_vals - min_vals
        )
    
    calcTime = 0
    abort = False
    sortiePhreeqCtotal = pd.DataFrame()
    for index in commMtrxPart.index:

        if abort : break
        if pcDict['maillesChargeGeom']: specieChargeMailleListe = speciesToNode(index, pcDict['maillesChargeGeom'], pcDict['speciesChargeGeometry']) # associe une liste d'espèces de contre-charge en fonction de la maille
        else: specieChargeMailleListe =  None
        
        try:
            ref = time.perf_counter()
            phreeqc.run_string(nodeSpeciation(commMtrxPart, commMtrx_primSpecies, pcDict['speciesCharge'],specieChargeMailleListe, pcDict['speciationCharge'],pcDict['solMod']))

            calcTime += time.perf_counter() - ref
            
        except Exception as e: # a changer
            
            if pcDict['speciationCharge']:
                start = 1
                if pcDict['speciesChargeGeometry']:
                    scriptWarnings +=f"PhreeqC : node n°{index}, time step n°{pcDict['lStep']+1}, t={pcDict['tStep']}{pcDict['timeUnit']}, PID={os.getpid()} : assuming {specieChargeMailleListe[0][0]} as counter-charge species :\n {e}\n"
                elif pcDict['speciesCharge']:
                    scriptWarnings +=f"PhreeqC : node n°{index}, time step n°{pcDict['lStep']+1}, t={pcDict['tStep']}{pcDict['timeUnit']}, PID={os.getpid()} : assuming {pcDict['speciesCharge'][0]} as counter-charge species :\n {e}\n"
            else:
                start = 0
                scriptWarnings +=f"PhreeqC : node n°{index}, time step n°{pcDict['lStep']+1}, t={pcDict['tStep']}{pcDict['timeUnit']}, PID={os.getpid()} : assuming no counter-charge species :\n {e}\n"

            if pcDict['speciesChargeGeometry']:
                for k,spc in enumerate(specieChargeMailleListe[0][start:], start=start): 
                    try:
                        warnings +=1
                        scriptWarnings += f"PhreeqC : node n°{index}, time step n°{pcDict['lStep']+1}, t={pcDict['tStep']}{pcDict['timeUnit']}, PID={os.getpid()} : testing {spc} as geometrical dependent counter-balance species ...\n"
                        ref = time.perf_counter()
                        phreeqc.run_string(nodeSpeciation(commMtrxPart, commMtrx_primSpecies, None, spc, True, pcDict['solMod']))
                        calcTime += time.perf_counter() - ref
                        scriptWarnings +=  f"PhreeqC : node n°{index}, time step n°{pcDict['lStep']+1}, t={pcDict['tStep']}{pcDict['timeUnit']}, PID={os.getpid()} : success !\n"
                        break
                    except Exception as r:
                        scriptWarnings += f"PhreeqC : node n°{index}, time step n°{pcDict['lStep']+1}, t={pcDict['tStep']}{pcDict['timeUnit']}, PID={os.getpid()} : {r}\n"
                    if spc == specieChargeMailleListe[0][-1]:
                        scriptWarnings += f"PhreeqC : node n°{index}, time step n°{pcDict['lStep']+1}, t={pcDict['tStep']}{pcDict['timeUnit']}, PID={os.getpid()} : Fatal PhreeqC error. Aborted PhreeqC batch :\n"
                        scriptWarnings += nodeSpeciation(commMtrxPart, commMtrx_primSpecies, None, spc, True, pcDict['solMod'])
                        print('Speciation batch aborted. See warning.log file.')
                        abort = True
                        
            elif pcDict['speciesCharge']:
                for k,spc in enumerate(pcDict['speciesCharge'][start:], start=start):
                    try:
                        warnings +=1
                        scriptWarnings +=f"PhreeqC : node n°{index}, time step n°{pcDict['lStep']+1}, t={pcDict['tStep']}{pcDict['timeUnit']}, PID={os.getpid()} : testing {spc} as counter-balance species ...\n"
                        ref = time.perf_counter()
                        phreeqc.run_string(nodeSpeciation(commMtrxPart, commMtrx_primSpecies, spc, None, True,pcDict['solMod']))
                        calcTime += time.perf_counter() - ref
                        scriptWarnings +=f"PhreeqC : node n°{index}, time step n°{pcDict['lStep']+1}, t={pcDict['tStep']}{pcDict['timeUnit']}, PID={os.getpid()} : success !\n"
                        break
                    except Exception as r: 
                        scriptWarnings +=f"PhreeqC : node n°{index}, time step n°{pcDict['lStep']+1}, t={pcDict['tStep']}{pcDict['timeUnit']}, PID={os.getpid()} : {r}\n"
                        
                    if spc == pcDict['speciesCharge'][-1]:
                        scriptWarnings +=f"PhreeqC : node n°{index}, time step n°{pcDict['lStep']+1}, t={pcDict['tStep']}{pcDict['timeUnit']}, PID={os.getpid()} : Fatal PhreeqC error. Aborted PhreeqC batch :\n"
                        scriptWarnings += nodeSpeciation(commMtrxPart, commMtrx_primSpecies, None, spc, True, pcDict['solMod'])
                        print('Speciation batch aborted. See warning.log file.')
                        abort = True
                                
            else:
                scriptWarnings +=f"PhreeqC : node n°{index}, time step n°{pcDict['lStep']+1}, t={pcDict['tStep']}{pcDict['timeUnit']}, PID={os.getpid()} : Fatal PhreeqC error. Aborted PhreeqC batch :\n"
                scriptWarnings += nodeSpeciation(commMtrxPart, commMtrx_primSpecies, pcDict['speciesCharge'][0],specieChargeMailleListe, pcDict['speciationCharge'],pcDict['solMod'])
                print('Speciation batch aborted. See warning.log file.')
                abort = True
        
        if phreeqc.get_selected_output_array() and not abort:
            
            out = pd.DataFrame([phreeqc.get_selected_output_array()[-1]],  columns=phreeqc.get_selected_output_array()[0])
            sortiePhreeqCtotal = pd.concat([sortiePhreeqCtotal, out], ignore_index=True) 
            # if pcDict["firstStepEquilibrium"]==False:
            # if index == 1 : print(sortiePhreeqCtotal)
            
            
        elif not phreeqc.get_selected_output_array():
            warnings+=1
            scriptWarnings +=f"PhreeqC : node n°{index}, time step n°{pcDict['lStep']+1}, t={pcDict['tStep']}{pcDict['timeUnit']}, PID={os.getpid()} : no PhreeqC ouput.\n"
            scriptWarnings += nodeSpeciation(commMtrxPart, commMtrx_primSpecies, None, None, pcDict['speciationCharge'], pcDict['solMod'])
            abort = True
    # print(sortiePhreeqCtotal)
    if abort:
        return pd.DataFrame(),pd.DataFrame(),pd.DataFrame(),pd.DataFrame(),warnings, scriptWarnings, abort,calcTime 
    else:

        sortiePhreeqCtotal.index = commMtrxPart.index
        
        colonnesSpeciation  = [] # headers with phreeqc formalism, in piccts order
        if pcDict['primarySpecies'][1]: colonnesSpeciation += [spc if spc in (pcDict['primarySpecies'][1]) else f"m_{spc}(mol/kgw)" for spc in systemSpeciation if spc not in ['pH','pe'] ]
        else: colonnesSpeciation += [f"m_{spc}(mol/kgw)" for spc in systemSpeciation if spc not in ['pH','pe']]
        if 'pH' in systemSpeciation and 'pe' in systemSpeciation : colonnesSpeciation += ['pH','pe']

        commMtrxSpct = sortiePhreeqCtotal[colonnesSpeciation].copy()
        commMtrxSpct.columns = systemSpeciation
        
        # coupled parameters
        # commMtrxSpct = pd.concat([commMtrxPart[pcDict['coord']],commMtrxSpct], axis = 1)
        # commMtrxSpct = commMtrxSpct[list(commMtrxPart.columns)]
        

        col = [c for c in sortiePhreeqCtotal.columns if c not in colonnesSpeciation]
        colFinale = col  +  colonnesSpeciation
        sortiePhreeqCtotal[colFinale]
        
        # sortiePhreeqCtotal = pd.concat([commMtrxPart[pcDict['coord']], sortiePhreeqCtotal], axis=1)
        
        if pcDict['acidicTrspt']:
            AcidicEcho = sortiePhreeqCtotal[['m_H+(mol/kgw)', 'm_OH-(mol/kgw)','pH']].copy()
            AcidicEcho.columns = ['H+','OH-','pH']
            AcidicEcho.index = commMtrxPart.index
            return commMtrxSpct, AcidicEcho, commMtrx_primSpecies, sortiePhreeqCtotal, warnings, scriptWarnings, abort,calcTime 
        
        else:
            return commMtrxSpct, pd.DataFrame(), commMtrx_primSpecies, sortiePhreeqCtotal, warnings, scriptWarnings, abort,calcTime 

    
def spct(centralDict):
    startPhreeqC = time.time()
    print("PhreeqC", end=" ", flush=True)
    
    initialColumns = list(centralDict['commMtrx'].columns)
    
    if centralDict['partialCrossDependencies'] and centralDict['partialCrossDependencies']['chem']:
        crossDepMtrx = pd.DataFrame()
        crossDepMtrx = centralDict['commMtrx'][centralDict['partialCrossDependencies']['chem']].copy()
        commMtrx = centralDict["commMtrx"].drop(columns=centralDict['partialCrossDependencies']['chem'])
    else: commMtrx = centralDict["commMtrx"].copy()

    
    commMtrxCoord = commMtrx[centralDict['coord']].copy()
    commMtrx = commMtrx.drop(columns=centralDict['coord'])
    
    chunk_size = int(np.ceil(len(commMtrx) / centralDict['PIDnbr']))
    # print(commMtrx)
    commMtrxSplit = [commMtrx.iloc[i:i + chunk_size] for i in range(0, len(commMtrx), chunk_size)]

    if centralDict['renouvellement'] and centralDict['renouvellement'][0] <= centralDict["tStep"] :
        rnvlmt = True
        del centralDict['renouvellement'][0]
    else: rnvlmt = False
    
    with concurrent.futures.ProcessPoolExecutor(max_workers=centralDict['PIDnbr']) as executor:
        futures = []
        for chunk in commMtrxSplit:
            futures.append(executor.submit(speciationPhreeqC, centralDict, chunk, rnvlmt))
            

    results = []
    results = [f.result() for f in futures]

    df1, df2, df3, df4, intgr, strg, Bool, intgr2 = zip(*results)

    commMtrxSpct = pd.concat(df1, ignore_index=True)
    AcidicEcho = pd.concat(df2, ignore_index=True)
    commMtrx_primSpecies = pd.concat(df3, ignore_index=True)
    sortiePhreeqCtotal = pd.concat(df4, ignore_index=True)
    totalWarnings = sum(intgr)
    prcsTime = sum(intgr2)
    totalScriptWarnings = "".join(strg)
    
    
    if totalWarnings:
        with open("warning.log", "a") as warningLog:
            warningLog.write(f"PhreeqC, time = {centralDict['tStep']}{centralDict['timeUnit']}, time step n°{centralDict['lStep']+1} : the {totalWarnings} following warnings occured ...\n")
            warningLog.write(f"{totalScriptWarnings}\n")
    if any(Bool):
        print('Fatal PhreeqC error. Aborting run.')
        sys.exit()
    
    if centralDict['partialCrossDependencies']['chem']:
        commMtrxSpct = pd.concat([commMtrxSpct,crossDepMtrx], axis=1)
    
    # print(commMtrxSpct)
    commMtrxSpct = pd.concat([commMtrxSpct,commMtrxCoord], axis=1)
    sortiePhreeqCtotal = pd.concat([commMtrxCoord,sortiePhreeqCtotal], axis=1)
    commMtrxSpct = commMtrxSpct[initialColumns]
    
    centralDict.update({
        "commMtrx": commMtrxSpct,
        "AcidicEcho": AcidicEcho,
        "PhreeqCClockTime": centralDict['PhreeqCClockTime']  + time.time() - startPhreeqC,
        "PhreeqCPrcsTime": centralDict['PhreeqCPrcsTime'] + prcsTime,
        })

    if centralDict["firstStepEquilibrium"]==True:
        commMtrx_primSpecies.to_csv(os.path.join(centralDict['paths']['primSpecies'], f"PrimarySpecies_{centralDict['lStep']}.txt"), index=False, header=True, sep='\t')
        sortiePhreeqCtotal.to_csv(os.path.join(centralDict['paths']['outputSpeciation'], f"PhreeqC_{centralDict['lStep']}.txt"), index=False, header=True, sep='\t')
    else:
        commMtrx_primSpecies.to_csv(os.path.join(centralDict['paths']['primSpecies'], f"PrimarySpecies_{centralDict['lStep']+1}.txt"), index=False, header=True, sep='\t')
        sortiePhreeqCtotal.to_csv(os.path.join(centralDict['paths']['outputSpeciation'], f"PhreeqC_{centralDict['lStep']+1}.txt"), index=False, header=True, sep='\t')

    print(f"({writeTime((time.time() - startPhreeqC))})") 

    return centralDict 