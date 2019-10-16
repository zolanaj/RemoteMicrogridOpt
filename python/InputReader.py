# -*- coding: utf-8 -*-
"""
This module contains a series of functions that allow for the user
to import input data from a series of csv files (either comma or
space delimited).  
"""

import csv
#import pandas
from DieselGenerator import DieselGenerator
from PVArray import PVArray
from Battery import Battery

  
def ReadTechnologies(filename,soc_start=0.8,num_periods=8760):
    """Reads an input file and returns a list of the technologies available for 
    purchase.
    filename -- input file name
    soc_start -- the starting state of charge for any batteries purchased
    num_periods -- lenght of time horizon, in periods.
    
    retval -- a list of DieselGenerator, PVArray, and Battery objects.
    """
    technologies =[]
    tech_file = open(filename, 'rU')
    tech_reader = csv.reader(tech_file)
    keys = next(tech_reader)
    keys = [key.strip() for key in keys]
    for line in tech_reader:
        attributes = {}
        #print ', '.join(line)
        attributes["name"] = line[0]
        for i in range(1,len(line)):
            attributes[keys[i]] = line[i]
        #Make the appropriate technology based on the technology type. The
        #index given in the first column, "tech", is what we need.  Because the
        #notation is a letter-number combination, e.g. "G1", we use the
        #first character as a guide for the technology we create.
        if attributes["name"][0] == 'G':
            #make a DieselGenerator
            technologies.append(DieselGenerator(attributes, 0))
            #technologies[-1].SetFuelCost(fuel_cost)
        elif attributes["name"][0] == 'S':
            #make a PVArray
            technologies.append(PVArray(attributes, 0))
        elif attributes["name"][0] == 'B':
            #make a Battery
            technologies.append(Battery(attributes, 0) )
            technologies[-1].SetSOC(soc_start)
        else: assert False, "Error: invalid technology type."
    tech_file.close()
    return technologies
    
def ReadMultipleLoads(filename, overage, time_horizon):
    """Returns a dictionary of the load by period,
    by scenario, given a flat file as input.  Assumes a single
    PV array as given by the inputs."""
    loads = {}
    load_file = open(filename, 'rU')
    load_reader = csv.reader(load_file)
    headers = next(load_reader)
    headers = [head.strip() for head in headers]
    for key in headers[1:]:
        loads[key] = []
    for line in load_reader:
        for idx in range(1,len(line)):
            if len(loads[headers[idx]]) < time_horizon:
                loads[headers[idx]].append(max(5000.0,float(line[idx])*(1.0+overage)))
    return loads, headers[1:]

#def ReadMultiplePVInputs(filename, time_horizon):
#    """Returns a dictionary of the PV availability by period,
#    by scenario, given a flat file as input.  Assumes a single
#    PV array as given by the inputs."""
#    pvs = {}
#    pv_file = open(filename, 'rU')
#    pv_reader = csv.reader(pv_file)
#    headers = next(pv_reader)
#    headers = [head.strip() for head in headers]
#    for key in headers[1:]:
#        pvs[key] = []
#    for line in pv_reader:
#        for idx in range(1,len(line)):
#            if len(pvs[headers[idx]]) < time_horizon:
#                pvs[headers[idx]].append(float(line[idx]))
#    return pvs

def ReadMultiplePVInputs(filename, time_horizon):
    """Returns a dictionary of the PV availability by scenario,
    by type, by period, given a flat file as input.  Allows for
    multiple PV Arrays.  Input file must follow the format
    LL*.S*.TT*   X where * is a wildcard, LL is the scenario,
    S is the PV system type, TT is a time index, and X is a 
    floating point number (in Watts)."""
    pvs = {}
    pv_file = open(filename, 'rU')
    pv_reader = csv.reader(pv_file)
    for line in pv_reader:
        splitline = line[0].split(" ")
        inputs = splitline[0].split(".")
        if inputs[0] not in pvs.keys():  #add LL
            pvs[inputs[0]] = {}
        if inputs[1] not in pvs[inputs[0]].keys():  #Add S
            pvs[inputs[0]][inputs[1]] = {}
        pvs[inputs[0]][inputs[1]][inputs[2]] = float(splitline[-1])  #Add TT and data
    return pvs
        
def ReadPVSystemSpecs(filename):
    """Returns a dictionary of the PV System specs by array type."""
    pv_specs = {}
    pv_file = open(filename, 'rU')
    pv_reader = csv.reader(pv_file)
    headers = next(pv_reader)
    for line in pv_reader:
        name = line[0]
        pv_specs[name] = {}
        pv_specs[name]["name"] = name
        for idx in range(1,len(headers)):
            try: pv_specs[name][headers[idx]] = float(line[idx])
            except ValueError: pv_specs[name][headers[idx]] = line[idx]
    return pv_specs   

 
def ReadMultipleFuelCosts(filename, time_horizon):
    """Returns a dictionary of the fuel cost by period,
    by scenario, given a flat file as input."""
    costs = {}
    cost_file = open(filename, 'rU')
    cost_reader = csv.reader(cost_file)
    headers = next(cost_reader)
    headers = [head.strip() for head in headers]
    for key in headers[1:]:
        costs[key] = []
    for line in cost_reader:
        for idx in range(1,len(line)):
            if len(costs[headers[idx]]) < time_horizon:
                costs[headers[idx]].append(float(line[idx]))
    return costs   
    
def ReadTechMaxes(filename):
    """Returns a dictionary of the maximum allowable purcahses
    by technology, by scenario, given a flat file as input."""
    maxes = {}
    max_file = open(filename, 'rU')
    max_reader = csv.reader(max_file)
    headers = next(max_reader)
    headers = [head.strip() for head in headers]
    for line in max_reader:
        key = line[0].strip()
        maxes[key] = {}
        for idx in range(1,len(line)):
            maxes[key][headers[idx]] = int(line[idx])
    return maxes 
    
def ReadLoadsAndFuelCosts(filename, overage, time_horizon):
    loads = []
    fuel_costs = []
    load_file = open(filename, 'rU')
    load_reader = csv.reader(load_file)
    next(load_reader)
    for line in load_reader:
        #time_index = line[0]
        load = line[1]
        loads.append(float(load)*(1.0+overage))
        fuel_cost = line[2]
        fuel_costs.append(float(fuel_cost))
    #Update load to account for overage
    load_file.close()
    return loads[:time_horizon], fuel_costs[:time_horizon]
    
def ReadLoadPV(filename, overage, time_horizon):
    loads = []
    pvs = []
    load_file = open(filename, 'rU')
    load_reader = csv.reader(load_file)
    next(load_reader)
    for line in load_reader:
        #time_index = line[0]
        time, load, pv = line
        loads.append(float(load)*(1.0+overage))
        pvs.append(float(pv))
    #Update load to account for overage
    load_file.close()
    return loads[:time_horizon], pvs[:time_horizon]   
    
def ReadLoadCostTemp(filename, overage, time_horizon):
    output = {}
    load_file = open(filename, 'rU')
    load_reader = csv.reader(load_file)
    headers = next(load_reader)
    headers = [head.strip() for head in headers]
    for item in headers:
        output[item] = []
    for line in load_reader:
        for idx, val in enumerate(line):
            if headers[idx] == "d_p":
                output[headers[idx]].append(float(val)*(1.0+overage))
            elif idx >= 1:
                output[headers[idx]].append(float(val))
    loads = output["d_p"][:time_horizon]
    fuel_costs = output["delta_fuel"][:time_horizon]
    if "temp_t" in output.keys():
        temps = output["temp_t"][:time_horizon]
    else: 
        temps = None
    if "solar" in output.keys():
        pvs = output["temp_t"][:time_horizon]
    else: 
        pvs = None
    load_file.close()
    #print temps
    return loads, fuel_costs, temps, pvs
    
def ReadScalars(filename):
    scalars = {}
    scalar_file = open(filename,'rU')
    scalar_reader = csv.reader(scalar_file)
    for line in scalar_reader:
        splitline = line[0].split()
        scalars[splitline[0]] = float(splitline[1])
    #print "Scalars", scalars
    scalar_file.close()
    return scalars

def ReadPVFile(pv_filename, technologies):
    pv_avails = {}
    pv_file = open(pv_filename, 'rU')
    pvs_in = csv.reader(pv_file)
    for line in pvs_in:
        splitline = line[0].split()
        tech = splitline[0][:2]
        if tech not in pv_avails.keys():
            pv_avails[tech] = [float(splitline[1])]
        else:
            pv_avails[tech].append(float(splitline[1]))
    for tech in technologies:
        if tech.GetType() == "PVArray":
            tech.SetPVAvailability(pv_avails[tech.GetName()])
    return technologies
    
def ReadTemps(temp_filename):
    temps = []
    temp_file = open(temp_filename, 'rU')
    temps_in = csv.reader(temp_file)
    for line in temps_in:
        splitline = line.split()  
        temps.append(float(splitline[1]))    
    return temps 
    
def ReadInitialSols1Index(filename):
    """Takes as input a GAMS-style csv input file with the notation A*,
    where * is a wildcard and A is any index, followed by two spaces and 
    a number.  This is read and returned as a dictionary with nested keys
    that is used by the CPLEX solver."""
    d = {}
    infile = open(filename, 'rU')
    reader = csv.reader(infile)
    for line in reader:
        splitline = line[0].split()
        d[splitline[0]] = float(splitline[-1])  #Add A and data
    return d    
    
def ReadInitialSols2Index(filename,style='i'):
    """Takes as input a GAMS-style csv input file with the notation A*.B*,
    where * is a wildcard and A and B are any index, followed by two spaces and 
    a number.  This is read and returned as a dictionary with nested keys
    that is used by the CPLEX solver."""
    d = {}
    infile = open(filename, 'rU')
    reader = csv.reader(infile)
    for line in reader:
        splitline = line[0].split()
        params = splitline[0].split(".")
        if params[0] not in d.keys(): d[params[0]] = {}  #Add A
        if style == 'i':
            d[params[0]][params[1]] = int(.5+float(splitline[-1])) #add C and data
        else:
            d[params[0]][params[1]] = float(splitline[-1]) #add C and data
    return d
        
def ReadInitialSols3Index(filename,style='i'):
    """Takes as input a GAMS-style csv input file with the notation A*.B*.C*,
    where * is a wildcard and A through C are any index, followed by two spaces 
    and a number.  This is read and returned as a dictionary with nested keys
    that is used by the CPLEX solver."""
    d = {}
    infile = open(filename, 'rU')
    reader = csv.reader(infile)
    for line in reader:
        splitline = line[0].split()
        params = splitline[0].split(".")
        if params[0] not in d.keys():  #add A
            d[params[0]] = {}
        if params[1] not in d[params[0]].keys():  #add B
            d[params[0]][params[1]] = {}
        if style == 'i':
            d[params[0]][params[1]][params[2]] = int(.5+float(splitline[-1])) #add C and data
        else:
            d[params[0]][params[1]][params[2]] = float(splitline[-1]) #add C and data
    return d
    
def ReadInitialSols4Index(filename,style='f'):
    """Takes as input a GAMS-style csv input file with the notation A*.B*.C*.D*,
    where * is a wildcard and A through D are any index, followed by two spaces 
    and a number.  This is read and returned as a dictionary with nested keys
    that is used by the CPLEX solver."""
    d = {}
    infile = open(filename, 'rU')
    reader = csv.reader(infile)
    for line in reader:
        splitline = line[0].split()
        params = splitline[0].split(".")
        if params[0] not in d.keys():  #add A
            d[params[0]] = {}
        if params[1] not in d[params[0]].keys():  #add B
            d[params[0]][params[1]] = {}
        if params[2] not in d[params[0]][params[1]].keys():  #add C
            d[params[0]][params[1]][params[2]] = {}
        if style == 'i':
            d[params[0]][params[1]][params[2]][params[3]] = int(.5+float(splitline[-1])) #add D and data
        else:
            d[params[0]][params[1]][params[2]][params[3]] = float(splitline[-1]) #add D and data
    return d
"""    
def GetInitialSolution
"""   
def GetSpaceDelimInputs(filename,overage=0.3):
    d = {}
    infile = open(filename, 'rU')
    reader = csv.reader(infile)
    headline = next(reader)
    header_split = headline[0].split()
    #print header_split
    for x in header_split:
        d[x] = []
    for line in reader:
        splitline = line[0].split()
        #print splitline
        for i in range(1,len(splitline)):
            if header_split[i+1] == "d_P":
                try: d[header_split[i+1]].append( (1+overage)*float(splitline[i]) )
                except ValueError: d[header_split[i+1]].append( (1+overage)*float(splitline[i][:-1]) )
            else:
                try: d[header_split[i+1]].append( float(splitline[i]) )
                except ValueError: d[header_split[i+1]].append( float(splitline[i][:-1]) )
    return d
    
def ReadAMPLTechMaxes(filename):
    d = {}
    infile = open(filename, 'rU')
    reader = csv.reader(infile)
    for line in reader:
        splitline = line[0].split()
        d[splitline[0]] = int(splitline[1])
    return d
    
def GetVarsList2(d,varname,scenario):
    vals = []
    for key1 in d[scenario].keys():
        vals.append((varname+key1,d[scenario][key1]))
    return vals
    
def GetVarsList3(d,varname,scenario):
    vals = []
    for key1 in d[scenario].keys():
        for key2 in d[scenario][key1].keys():
            vals.append((varname+key1+"."+key2,d[scenario][key1][key2]))
    return vals
    
def GetVarsList4(d,varname,scenario):
    vals = []
    for key1 in d[scenario].keys():
        for key2 in d[scenario][key1].keys():
            for key3 in d[scenario][key1][key2].keys():
                vals.append((varname+key1+"."+key2+"."+key3,d[scenario][key1][key2][key3]))
    return vals
    
def GetIFSVals(ifs_dir,scenario,binit=0.8):
    W_file = ifs_dir+"WW.csv"
    X_file = ifs_dir+"X.csv"
    L_file = ifs_dir+"L.csv"
    I_plus_file = ifs_dir+"I_plus.csv"
    I_minus_file = ifs_dir+"I_minus.csv"
    P_in_file = ifs_dir+"P_in.csv"
    P_out_file = ifs_dir+"P_out.csv"
    B_minus_file = ifs_dir+"B_minus.csv"
    B_plus_file = ifs_dir+"B_plus.csv"
    Z_minus_file = ifs_dir+"Z_minus.csv"
    Z_plus_file = ifs_dir+"Z_plus.csv"
    P_solar_file = ifs_dir+"P_solar.csv"
    B_soc_file = ifs_dir+"B_soc.csv"
    F_tilde_file = ifs_dir+"F_tilde.csv" 
    G_file = ifs_dir+"G.csv"
    ifs = {}
    ifs["W"] = GetVarsList3(ReadInitialSols3Index(W_file,style='i'),"W_",scenario)
    ifs["X"] = GetVarsList2(ReadInitialSols2Index(X_file,style='i'),"X_",scenario)
    ifs["L"] = GetVarsList3(ReadInitialSols3Index(L_file,style='f'),"L_",scenario)
    ifs["I_plus"] = GetVarsList4(ReadInitialSols4Index(I_plus_file,style='f'),"I_plus_",scenario)
    ifs["I_minus"] = GetVarsList4(ReadInitialSols4Index(I_minus_file,style='f'),"I_minus_",scenario)
    ifs["Y_plus"] = GetVarsList4(ReadInitialSols4Index(I_plus_file,style='f'),"Y_plus_",scenario)
    ifs["Y_minus"] = GetVarsList4(ReadInitialSols4Index(I_minus_file,style='f'),"Y_minus_",scenario)
    ifs["P_plus"] = GetVarsList4(ReadInitialSols4Index(P_in_file,style='f'),"P_plus_",scenario)
    ifs["P_minus"] = GetVarsList4(ReadInitialSols4Index(P_out_file,style='f'),"P_minus_",scenario)
    ifs["B_minus"] = GetVarsList4(ReadInitialSols4Index(B_minus_file,style='i'),"B_minus_",scenario)
    ifs["B_plus"] = GetVarsList4(ReadInitialSols4Index(B_plus_file,style='i'),"B_plus_",scenario)
    ifs["Z_minus"] = GetVarsList4(ReadInitialSols4Index(Z_minus_file,style='f'),"Z_minus_",scenario)
    ifs["Z_plus"] = GetVarsList4(ReadInitialSols4Index(Z_plus_file,style='f'),"Z_plus_",scenario)
    ifs["P_PV"] = GetVarsList3(ReadInitialSols3Index(P_solar_file,style='f'),"P_PV_",scenario)
    ifs["B_soc"] = GetVarsList4(ReadInitialSols4Index(B_soc_file,style='f'),"B_soc_",scenario)
    ifs["F_tilde"] = GetVarsList2(ReadInitialSols2Index(F_tilde_file,style='f'),"F_tilde_",scenario)
    ifs["G"] = GetVarsList4(ReadInitialSols4Index(G_file,style='i'),"G_",scenario)
    ifs["socstart"] = [] 
    for ind, val in ifs["W"]:
        if "B" in ind and val == 1:
            ifs["socstart"].append(("soc_start_"+ind[2:],binit))
    return ifs
        
    
    
def GetInputs(scalar_filename,load_filename,pv_filename,pv_spec_filename,
            fc_filename,tech_filename,max_filename,cut_filename,time_horizon,
            mincap_filename = None, W_filename = None, X_filename = None):
    scalars = ReadScalars(scalar_filename)
    #if 'nu' in scalars.keys(): nu = float(scalars['nu'])        #nu - ratio of operation time to model time
    #if 'delfuel' in scalars.keys(): delfuel = float(scalars['delfuel'])
    overage = float(scalars['k_bar']) 
    soc_start = float(scalars['b_init'])
    tech_maxes = ReadTechMaxes(max_filename)
    technologies = ReadTechnologies(tech_filename, 
                                soc_start,time_horizon)
    loads, scenarios = ReadMultipleLoads(load_filename,overage,time_horizon)
    fuel_costs = ReadMultipleFuelCosts(fc_filename,time_horizon)
    pv_avail = ReadMultiplePVInputs(pv_filename, time_horizon)
    pv_specs = ReadPVSystemSpecs(pv_spec_filename)
    cuts = ReadTechMaxes(cut_filename)
    if mincap_filename != None: mincaps = ReadScalars(mincap_filename)
    else: mincaps = None
    if W_filename != None: W_init = ReadInitialSols3Index(W_filename)
    else: W_init = None
    if X_filename != None: X_init = ReadInitialSols2Index(X_filename)
    else: X_init = None
    #removing batteries we don't want - equivalent of subset given in GAMS                                                            
    for pv in pv_specs.keys():
        if pv_specs[pv]["include"] > 0.5:
            technologies.append( PVArray(
                    pv_specs[pv],
                    0
                    )
                )
            technologies[-1].ConvertPVAvailability(len(loads))
    #give initial parameters to batteries, which aren't given in 
    #the technology inputs file.
    soc_min = float(scalars['soc_min'])  #min soc of batteries
    soc_max = float(scalars['soc_max'])  #max soc of batteries
    i_min = max(0.0,float(scalars['min_current']))  #minimum current for batteries                                                                                                            
    for tech in technologies:
        if tech.GetType() == "Battery":
            tech.SetMinSOC(soc_min)
            tech.SetMaxSOC(soc_max)
            tech.SetMinCurrent(i_min)
            #tech.SetTechMax(2)  #placeholder for max # of batteries
        #elif tech.GetType() == "PVArray":
        #    tech.SetTechMax(int(scalars["solar_max"]))
    #set maximum number of technologies, by PV Array.
    return {"technologies":technologies, 
            "scalars":scalars, "loads":loads, "pv_avail":pv_avail, 
            "fuel_costs":fuel_costs, "tech_maxes":tech_maxes, 
            "scenarios":scenarios, "cuts":cuts, 
            "mincaps": mincaps,
             "X_init": X_init,
             "W_init": W_init}

def setBounds(inputs,scenario,time_horizon=8760):
    for tech in inputs["technologies"]:
        if tech.GetType() != "PVArray":
            tech.SetTechMax(inputs["tech_maxes"][scenario][tech.GetName()])
        else: 
            if scenario[-2] == "l" or int(scenario[-2:]) <= 14:  tech.SetPVAvailability(inputs["pv_avail"][scenario][tech.GetName()])
            else: tech.SetPVAvailability(inputs["pv_avail"]["ll12"][tech.GetName()])
            tech.ConvertPVAvailability(time_horizon)
             
def GetAMPLInputs(scalar_filename,ampl_filename,pv_spec_filename,
            tech_filename,max_filename,cut_filename,time_horizon,
            mincap_filename = None, W_filename = None, X_filename = None):
    scalars = ReadScalars(scalar_filename)
    overage = float(scalars['k_bar']) 
    soc_start = float(scalars['b_init'])
    tech_maxes = GetAMPLTechMaxes(max_filename)
    technologies = ReadTechnologies(tech_filename, 
                                soc_start,time_horizon)
    ampl_in = GetSpaceDelimInputs(ampl_filename,overage)
    loads = ampl_in["d_P"]
    fuel_costs = ampl_in["delta_f"]
    pvs = ampl_in["gamma_t"]
    pv_specs = ReadPVSystemSpecs(pv_spec_filename)
    cuts = ReadTechMaxes(cut_filename)
    if mincap_filename != None: mincaps = ReadScalars(mincap_filename)
    else: mincaps = None
    if W_filename != None: W_init = ReadInitialSols3Index(W_filename)
    else: W_init = None
    if X_filename != None: X_init = ReadInitialSols2Index(X_filename)
    else: X_init = None
    #removing batteries we don't want - equivalent of subset given in GAMS                                                            
    for pv in pv_specs.keys():
        if pv_specs[pv]["include"] > 0.5:
            technologies.append( PVArray(
                    pv_specs[pv],
                    0
                    )
                )
            technologies[-1].SetPVAvailability(pvs)
    #give initial parameters to batteries, which aren't given in 
    #the technology inputs file.
    soc_min = float(scalars['soc_min'])  #min soc of batteries
    soc_max = float(scalars['soc_max'])  #max soc of batteries
    i_min = max(0.0,float(scalars['min_current']))  #minimum current for batteries                                                                                                            
    for tech in technologies:
        if tech.GetType() == "Battery":
            tech.SetMinSOC(soc_min)
            tech.SetMaxSOC(soc_max)
            tech.SetMinCurrent(i_min)
            #tech.SetTechMax(2)  #placeholder for max # of batteries
        elif tech.GetType() == "PVArray":
            tech.SetTechMax(int(scalars["solar_max"]))
    return {"technologies": technologies, 
             "scalars": scalars, 
             "loads": loads, 
             "pvs": pvs, 
             "fuel_costs": fuel_costs, 
             "tech_maxes": tech_maxes, 
             #"scenarios": scenarios, 
             "cuts": cuts,
             "mincaps": mincaps,
             "X_init": X_init,
             "W_init": W_init
             }  
             
#if __name__ == "__main__":
#    tm = GetAMPLTechMaxes("x_max.dat")
#    print tm
#    sc = ReadScalars("x_max.dat")
#    print sc