"""
Written by Alex Zolan 
Univariate Paritioning Model (U)

This module contains a collection of methods used to build and solve instances
of a Remote Microgrid Optimization Model, which we refer to as Model (U) in the
paper.  This model uses McCormick sub-envelopes using a univariate partitioning
scheme to approximate bilinear terms for the product of a battery's current and
its state of charge.

Note that we use a separate module, InputReader, for the 
purpose of reading inputs.  
"""
#import DieselGenerator, PVArray, Battery
import multiprocessing
import scipy
import InputReader 
import cplex
import time

#update working directory to current directory.
import os

#setting one global variable
local = True

def FractionOfSolar(tech,start_period,end_period):
    return scipy.sum([tech.CalculatePowerAvailable(t) 
            for t in range(start_period-1,end_period)]) / scipy.sum([
            tech.CalculatePowerAvailable(t) 
            for t in range(len(tech.GetPVAvailability())) ])

def FractionOfDemand(loads,start_period,end_period):
    return sum(loads[start_period-1:end_period]) / sum(loads)
    
def SolvePartition(inputs, scenario, start_period, end_period, 
        inv_mult=0.0, upper_bound = False, mipgap=0.05, 
        reslim=3600, test=False, design_mults = [0,0],
        IFS = None, batPartsMinus = [0.0,0.25,0.5,0.75,1.0], 
        batPartsPlus = [0.0,0.25,0.5,0.75,1.0], boundary_soc = 0.5,
        output=False, amortize_on_demand = True, min_PV=0, procs=1, mincap=0):
    """creates a full partition of the model.
    start_period -- starting period for the model
    end_period -- ending period for the model
    inv_mult -- price of SOC constraint violation (start = end = reset point)
    upper_bound -- indicator of whether the upper bound or lower
    bound is to be solved
    rhsmin -- input from the scenario - this is not yet in 
        InputReader
    batParts -- fraction of max current over which the McCormick underestimators
        for current operate.  Default values are evenly spaced on the [0,1] 
        interval but our input comes from the linear model solution and is the
        25th,50th, and 75th percentiles of current (as well as the min and max 
        allowable).
    
    retval -- cplex solution and objective value
    """
    numBatParts = len(batPartsPlus)-1
    technologies = inputs["technologies"]
    generators = [tech for tech in technologies if tech.GetType() == "DieselGenerator"]
    batteries = [tech for tech in technologies if tech.GetType() == "Battery"]
    pv_arrays = [tech for tech in technologies if tech.GetType() == "PVArray"]
    tech_maxes = inputs["tech_maxes"][scenario]
    scalars = inputs["scalars"]
    loads = scipy.array(inputs["loads"][scenario][
            start_period-1:end_period])
    pvs = inputs["pv_avail"][scenario]
    fuel_costs = scipy.array(inputs["fuel_costs"][scenario])
    tech_maxes = inputs["tech_maxes"][scenario]
    mingens = inputs["cuts"][scenario]['rhsmin']
    num_periods = end_period+1-start_period
    prob = 1.0*num_periods/8760
    for tech in technologies:
        if tech.GetType() != "PVArray":
            tech.SetTechMax(tech_maxes[tech.GetName()])
        else: 
            tech.SetPVAvailability(pvs[tech.GetName()])
            tech.ConvertPVAvailability(len(inputs["loads"][scenario]))
    W_init = inputs["W_init"]
    X_init = inputs["X_init"]    
    #print technologies[0].GetType()
    #print tech_maxes
    #initialize problem and suppress console output
    p = cplex.Cplex()
    p.set_error_stream(None)
    p.set_warning_stream(None)
    if not output:
        p.set_log_stream(None)
        p.set_results_stream(None)
    else: 
        p.set_log_stream(scenario+"-"+str(start_period)+"_log.txt")
        if upper_bound:
            resfilename = scenario+"-"+str(start_period)+"_UB.txt"
        else: 
            resfilename = scenario+"-"+str(start_period)+"_LB.txt"
        p.set_results_stream(resfilename)
        
    p.parameters.timelimit.set(reslim) 
    p.parameters.workmem.set(50)
    p.parameters.mip.tolerances.mipgap.set(mipgap)
    p.parameters.threads.set(procs)
    p.parameters.lpmethod.set(4)
    p.parameters.mip.strategy.rinsheur.set(10)
    p.parameters.mip.strategy.heuristicfreq.set(100)
    p.parameters.mip.strategy.file.set(2)
    #p.parameters.simplex.tolerances.feasibility.set(0.0001)
    p.parameters.emphasis.mip.set(2)
    p.parameters.mip.strategy.branch.set(1)
    def flatten(l):
        out = []
        for item in l:
            if isinstance(item, (list, tuple)):
                out.extend(flatten(item))
            else:
                out.append(item)
        return out
        
    def getName(name, *args):
        return '%s_%s'%(name, args)
        
    #define decision variable names
    def W_var(tech,k): return "W_"+tech.GetName()+".K"+str(k)
    def X_var(pv): return "X_"+pv.GetName()
    def P_minus_var(tech,k,t): return "P_minus_" + tech.GetName()+".K" + str(k) + ".TT" + str(t) 
    def P_plus_var(tech,k,t): return "P_plus_" + tech.GetName()+".K" + str(k) + ".TT" + str(t)
    def P_PV_var(pv,t): return "P_PV_"+pv.GetName()+".TT"+str(t)
    def F_tilde_var(t): return "F_tilde_TT"+str(t)
    def soc_start_var(bat,k): return "soc_start_"+bat.GetName() + ".K" + str(k)
    def B_soc_var(bat,k,t): 
        if t == start_period-1: 
            return soc_start_var(bat,k)
        return "B_soc_" + bat.GetName() + ".K" + str(k) + ".TT" + str(t) 
    def soc_end_var(bat,k): return B_soc_var(bat,k,end_period)
    def I_minus_var(bat,k,t): return "I_minus_" + bat.GetName() + ".K" + str(k) + ".TT" + str(t)
    def I_plus_var(bat,k,t): return "I_plus_" + bat.GetName() + ".K" + str(k) + ".TT" + str(t)
    def V_soc_var(bat,k,t): return "V_soc_" + bat.GetName() + ".K" + str(k) + ".TT" + str(t) 
    def Z_minus_var(bat,k,t): return "Z_minus_" + bat.GetName() + ".K" + str(k) + ".TT" + str(t)
    def Z_plus_var(bat,k,t): return "Z_plus_" + bat.GetName() + ".K" + str(k) + ".TT" + str(t)
    def B_minus_var(bat,k,t): return "B_minus_" + bat.GetName() + ".K" + str(k) + ".TT" + str(t)
    def B_plus_var(bat,k,t): return "B_plus_" + bat.GetName() + ".K" + str(k) + ".TT" + str(t)
    def G_var(gen,k,t): return "G_" + gen.GetName()+".K" + str(k) + ".TT" + str(t)
    def lambda_minus_var(bat,k,n,t): return "lambda_minus_"+bat.GetName() + ".K" + str(k) + ".N"+ str(n) + ".TT" + str(t)
    def lambda_plus_var(bat,k,n,t): return "lambda_plus_"+bat.GetName() + ".K" + str(k) + ".N"+ str(n) + ".TT" + str(t)
      
    W_names = flatten([["W_" + tech.GetName() + ".K" + str(i+1) 
            for i in range(tech_maxes[tech.GetName()])] 
            for tech in technologies 
                if tech.GetType() != "PVArray"])
    X_names = ["X_" + tech.GetName() for tech in technologies 
            if tech.GetType() == "PVArray"]
    P_minus_names = [[["P_minus_" + tech.GetName()+".K" + str(i+1) 
            + ".TT" + str(t) for i in range(tech_maxes[tech.GetName()])] 
            for tech in technologies 
                    if tech.GetType() != "PVArray"] 
            for t in range(start_period, end_period+1)]
    P_minus_names_gen = [[["P_minus_" + tech.GetName()+".K" + str(i+1) 
            + ".TT" + str(t) for i in range(tech_maxes[tech.GetName()])] 
            for tech in technologies 
                    if tech.GetType() == "DieselGenerator"] 
            for t in range(start_period, end_period+1)]
    P_plus_names = [[["P_plus_" + tech.GetName() + ".K" + str(i+1) 
             + ".TT" + str(t) for i in range(tech_maxes[tech.GetName()])] 
            for tech in technologies 
                    if tech.GetType() == "Battery"] 
            for t in range(start_period, end_period+1)]
    P_PV_names = [["P_PV_"+tech.GetName()+".TT"+str(t) 
            for tech in technologies 
                    if tech.GetType() == "PVArray"] 
            for t in range(start_period, end_period+1)]
    F_tilde_names = ["F_tilde_TT"+str(t) for t in range(
            start_period, end_period+1)]
    B_soc_names = [[["B_soc_" + tech.GetName() + ".K" + str(i+1) 
             + ".TT" + str(t) for i in range(tech_maxes[tech.GetName()])] 
            for tech in technologies 
                    if tech.GetType() == "Battery"] 
            for t in range(start_period, end_period+1)]
    I_minus_names = [[["I_minus_" + tech.GetName() + ".K" + str(i+1) 
             + ".TT" + str(t) for i in range(tech_maxes[tech.GetName()])] 
            for tech in technologies 
                    if tech.GetType() == "Battery"] 
            for t in range(start_period, end_period+1)]
    I_plus_names = [[["I_plus_" + tech.GetName() + ".K" + str(i+1) 
             + ".TT" + str(t) for i in range(tech_maxes[tech.GetName()])] 
            for tech in technologies 
                    if tech.GetType() == "Battery"] 
            for t in range(start_period, end_period+1)]
    Z_minus_names = [[["Z_minus_" + tech.GetName() + ".K" + str(i+1) 
             + ".TT" + str(t) for i in range(tech_maxes[tech.GetName()])] 
            for tech in technologies 
                    if tech.GetType() == "Battery"] 
            for t in range(start_period, end_period+1)]
    Z_plus_names = [[["Z_plus_" + tech.GetName() + ".K" + str(i+1) 
             + ".TT" + str(t) for i in range(tech_maxes[tech.GetName()])] 
            for tech in technologies 
                    if tech.GetType() == "Battery"] 
            for t in range(start_period, end_period+1)]
    B_minus_names = [[["B_minus_" + tech.GetName() + ".K" + str(i+1) 
             + ".TT" + str(t) for i in range(tech_maxes[tech.GetName()])] 
            for tech in technologies 
                    if tech.GetType() == "Battery"] 
            for t in range(start_period, end_period+1)]
    B_plus_names = [[["B_plus_" + tech.GetName() + ".K" + str(i+1) 
             + ".TT" + str(t) for i in range(tech_maxes[tech.GetName()])] 
            for tech in technologies 
                    if tech.GetType() == "Battery"] 
            for t in range(start_period, end_period+1)]    
    G_names = [[["G_" + tech.GetName() + ".K" + str(i+1) 
             + ".TT" + str(t) for i in range(tech_maxes[tech.GetName()])] 
            for tech in technologies 
                    if tech.GetType() == "DieselGenerator"] 
            for t in range(start_period, end_period+1)] 
    lambda_plus_names = [[[["lambda_plus_" + tech.GetName() + ".K" + str(i+1) 
            + ".N" + str(n+1) + ".TT" + str(t) 
            for i in range(tech_maxes[tech.GetName()])] 
            for tech in technologies  if tech.GetType() == "Battery"] 
            for n in range(numBatParts)]
            for t in range(start_period, end_period+1)] 
    lambda_minus_names = [[[["lambda_minus_" + tech.GetName() + ".K" + str(i+1) 
            + ".N" + str(n+1) + ".TT" + str(t) 
            for i in range(tech_maxes[tech.GetName()])] 
            for tech in technologies  if tech.GetType() == "Battery"] 
            for n in range(numBatParts)]
            for t in range(start_period, end_period+1)] 
    soc_start_names = [[soc_start_var(tech,(i+1))
            for i in range(tech_maxes[tech.GetName()])] 
            for tech in technologies 
                    if tech.GetType() == "Battery"] 
    soc_end_names = [[soc_end_var(tech,(i+1))
            for i in range(tech_maxes[tech.GetName()])] 
            for tech in technologies 
                    if tech.GetType() == "Battery"] 
    Jhat_flat = flatten(P_minus_names)
    Bhat_flat = flatten(P_plus_names) #batteries
    S_flat = flatten(P_PV_names)
    Ghat_flat = flatten(G_names)
    lamb_flat = flatten(lambda_plus_names)
    soc_start_flat = flatten(soc_start_names)
    num_Jhat = len(Jhat_flat)
    num_Bhat = len(Bhat_flat)
    num_S = len(S_flat)
    num_Ghat = len(Ghat_flat)
    num_lamb = len(lamb_flat)
    num_soc = len(soc_start_flat)
                
    #define parameters that will be used in calculating cost
    if amortize_on_demand:
        purchase_costs = ( FractionOfDemand(inputs["loads"][scenario],
                start_period,end_period) *
                scipy.array( flatten( [[tech.GetPurchaseCost() 
                for i in range(tech_maxes[tech.GetName()])] 
                for tech in technologies 
                if tech.GetType() != "PVArray"]) ) ) + design_mults[0]
        solar_costs = ( scipy.array([
                FractionOfSolar(tech,start_period,end_period) *
                tech.GetPurchaseCost() 
                for tech in technologies 
                if tech.GetType() == "PVArray"]) ) + design_mults[1]
    else:
        purchase_costs = prob*scipy.array(flatten([[tech.GetPurchaseCost() 
                for i in range(tech_maxes[tech.GetName()])] 
                for tech in technologies 
                    if tech.GetType() != "PVArray"])) + design_mults[0]    
        solar_costs = prob*scipy.array([tech.GetPurchaseCost() 
                for tech in technologies 
                if tech.GetType() == "PVArray"]) + design_mults[1]
    lambda_multipliers = scipy.array(flatten([[tech.GetMaxPower() 
            for i in range(tech_maxes[tech.GetName()])] 
            for tech in technologies 
                if tech.GetType() == "Battery"]))
    max_lambda_mult = lambda_multipliers.max()
    lambda_multipliers /= max_lambda_mult
    
    #lifecycle costs for I, Z variables (identical for charge and discharge)
    I_costs = [[[float(scalars["tau"])*tech.GetLifeCycleCost()*tech.GetASOC()/(2.*tech.GetReferenceCapacity())
            for i in range(tech_maxes[tech.GetName()])] 
            for tech in technologies 
                    if tech.GetType() == "Battery"] 
            for t in range(start_period, end_period+1)]
        
    Z_costs = [[[float(scalars["tau"])*tech.GetLifeCycleCost()*tech.GetDSOC()/(-2.*tech.GetReferenceCapacity())
            for i in range(tech_maxes[tech.GetName()])] 
            for tech in technologies 
                    if tech.GetType() == "Battery"] 
            for t in range(start_period, end_period+1)]
                    
    G_costs = [[[tech.GetLifeCycleCost() for i in range(tech_maxes[tech.GetName()])] 
            for tech in technologies 
                    if tech.GetType() == "DieselGenerator"] 
            for t in range(start_period, end_period+1)] 
        
    #creation of variables
    #W_bk -- purchase variables
    p.variables.add(obj=purchase_costs,
                    lb=scipy.zeros_like(purchase_costs), 
                    ub=scipy.ones_like(purchase_costs), 
                    types=[p.variables.type.binary]*len(purchase_costs), 
                    names= flatten(W_names) )
    #X_s -- Solar Purchase Costs
    p.variables.add(obj=solar_costs,
                    lb=scipy.zeros_like(solar_costs,dtype=float),
                    ub=[cplex.infinity]*len(solar_costs),  
                    types=[p.variables.type.continuous]*len(solar_costs), 
                    names= X_names )
    #P_minus_jkt -- Power out by variable by time period
    p.variables.add(obj=scipy.zeros(num_Jhat,dtype=float),
                    lb=scipy.zeros(num_Jhat,dtype=float),
                    ub=[cplex.infinity]*num_Jhat,
                    types=[p.variables.type.continuous]*num_Jhat,
                    names=flatten(P_minus_names) )
    #P_plus_jkt -- Power out by technology (no solar) by time period
    p.variables.add(obj=scipy.zeros(num_Bhat,dtype=float),
                    lb=scipy.zeros(num_Bhat,dtype=float),
                    ub=[cplex.infinity]*num_Bhat,
                    types=[p.variables.type.continuous]*num_Bhat,
                    names=flatten(P_plus_names) )
    #P_PV_st -- Power out by variable by time period
    p.variables.add(obj=scipy.zeros_like(S_flat,dtype=float),
                    lb=scipy.zeros(num_S,dtype=float),
                    ub=[cplex.infinity]*len(S_flat),
                    types=[p.variables.type.continuous]*num_S,
                    names=flatten(P_PV_names) )
    #F_tilde_t -- Fuel consumed (gal) by time period
    p.variables.add(obj=scipy.array(fuel_costs[start_period-1:end_period]),
                    lb=scipy.zeros_like(F_tilde_names,dtype=float),
                    ub=[cplex.infinity]*len(F_tilde_names),
                    types=[p.variables.type.continuous]*len(F_tilde_names),
                    names=F_tilde_names )
    #B_soc_bkt -- State of charge by battery by time period
    p.variables.add(obj=scipy.zeros_like(Bhat_flat,dtype=float),
                    lb=scipy.zeros(num_Bhat,dtype=float),
                    ub=scipy.ones(num_Bhat,dtype=float),
                    types=[p.variables.type.continuous]*num_Bhat,
                    names=flatten(B_soc_names) )
    #I_minus_bkt -- Current out by battery by time period
    p.variables.add(obj=scipy.array(flatten(I_costs)),
                    lb=scipy.zeros(num_Bhat,dtype=float),
                    ub=[cplex.infinity]*num_Bhat,
                    types=[p.variables.type.continuous]*num_Bhat,
                    names=flatten(I_minus_names) )
    #I_plus_bkt -- Current in by battery by time period
    p.variables.add(obj=scipy.array(flatten(I_costs)),
                    lb=scipy.zeros(num_Bhat,dtype=float),
                    ub=[cplex.infinity]*num_Bhat,
                    types=[p.variables.type.continuous]*num_Bhat,
                    names=flatten(I_plus_names) )
    #Z_minus_bkt -- Current out*State of charge by battery by time period
    p.variables.add(obj=scipy.array(flatten(Z_costs)),
                    lb=scipy.zeros(num_Bhat,dtype=float),
                    ub=[cplex.infinity]*num_Bhat,
                    types=[p.variables.type.continuous]*num_Bhat,
                    names=flatten(Z_minus_names) )
    #Z_plus_bkt -- Current in*State of charge by battery by time period
    p.variables.add(obj=scipy.array(flatten(Z_costs)),
                    lb=scipy.zeros(num_Bhat,dtype=float),
                    ub=[cplex.infinity]*num_Bhat,
                    types=[p.variables.type.continuous]*num_Bhat,
                    names=flatten(Z_plus_names) )
    #B_minus_bkt -- 1 if battery discharges, 0 o.w. by time period
    p.variables.add(obj=scipy.zeros_like(Bhat_flat,dtype=float),
                    lb=scipy.zeros(num_Bhat,dtype=float),
                    ub=scipy.ones(num_Bhat,dtype=float),
                    types=[p.variables.type.binary]*num_Bhat,
                    names=flatten(B_minus_names) )
    #B_plus_bkt -- 1 if battery charges, 0 o.w. by time period
    p.variables.add(obj=scipy.zeros_like(Bhat_flat,dtype=float),
                    lb=scipy.zeros(num_Bhat,dtype=float),
                    ub=scipy.ones(num_Bhat,dtype=float),
                    types=[p.variables.type.binary]*num_Bhat,
                    names=flatten(B_plus_names) )
    #G_gkt -- 1 if generator on, 0 o.w. by time period
    p.variables.add(obj=scipy.array(flatten(G_costs)),
                    lb=scipy.zeros(num_Ghat,dtype=float),
                    ub=scipy.ones(num_Ghat,dtype=float),
                    types=[p.variables.type.binary]*num_Ghat,
                    names=flatten(G_names) )
    #lambda_plus - 1 if SOC is in given in range and battery is charging, 0 o.w.
    p.variables.add(obj=scipy.zeros_like(lamb_flat,dtype=float),
                    lb=scipy.zeros(num_lamb,dtype=float),
                    ub=scipy.ones(num_lamb,dtype=float),
                    types=[p.variables.type.binary]*num_lamb,
                    names=flatten(lambda_plus_names) )
    #lambda_minus - 1 if SOC is in given in range and battery is charging, 0 o.w.
    p.variables.add(obj=scipy.zeros_like(lamb_flat,dtype=float),
                    lb=scipy.zeros(num_lamb,dtype=float),
                    ub=scipy.ones(num_lamb,dtype=float),
                    types=[p.variables.type.binary]*num_lamb,
                    names=flatten(lambda_minus_names) )
                                  
    #soc_start by scenario
    p.variables.add(obj=scipy.zeros(num_soc),
                    lb=[0.0]*num_soc,
                    ub=[1.0]*num_soc,
                    types=[p.variables.type.continuous]*num_soc,
                    names=flatten(soc_start_names) )
                    
    #inventory
    p.variables.add(obj=[inv_mult],
                    lb=[0.0],
                    ub=[cplex.infinity],
                    types=[p.variables.type.continuous],
                    names=['battery_inv'] )
    
    #prioritize branching
    orders = []
    for gen in generators:
        for k in range(tech_maxes[gen.GetName()]):
            orders.append((W_var(gen,k+1),1+int(gen.GetName()[1:]),p.order.branch_direction.up))
    for bat in batteries:
        for k in range(tech_maxes[bat.GetName()]):
            orders.append((W_var(bat,k+1),1,p.order.branch_direction.up)) 
    p.order.set(orders)
    #minimize cost
    p.objective.set_sense( p.objective.sense.minimize )

    #constraint 15a - power flow
    val = flatten([[tech.GetEfficiencyOut() 
            for i in range(tech_maxes[tech.GetName()])] 
            for tech in technologies 
            if tech.GetType() != "PVArray"])
    val.extend( [1.0 
            for tech in technologies 
            if tech.GetType() == "PVArray"]  )
    val.extend(flatten([[ -1.0
            for i in range(tech_maxes[tech.GetName()])] 
            for tech in technologies 
                if tech.GetType() == "Battery"]))
    p.linear_constraints.add(lin_expr=[
            cplex.SparsePair(flatten(P_minus_names[t]+
                    P_PV_names[t]+P_plus_names[t]), val)
            for t in range(num_periods)], 
            senses=['G' for t in range(num_periods)], 
            rhs=loads, 
            names=['Sys_Op_15a_'+str(t+start_period)
                for t in range(num_periods)])
    
    #constraint 15b - spinning reserve
    #eta_out*p_max*Bsoc for all batteries in time t
    val = flatten([[ tech.GetEfficiencyOut()*tech.GetMaxPower()*
            tech.GetSOC()
            for i in range(tech_maxes[tech.GetName()])] 
            for tech in technologies 
                if tech.GetType() == "Battery"])
    #p_max*G for all generators in time t
    val.extend(flatten([[ tech.GetMaxPower()
            for i in range(tech_maxes[tech.GetName()])] 
            for tech in technologies 
                if tech.GetType() == "DieselGenerator"]))
    #-P_minus for all generators in time t
    val.extend(scipy.ones(int(num_Ghat/num_periods)))
    #-k_sol*P_solar
    val.extend(-1.0*float(scalars["k_sol"])*scipy.ones(int(num_S/num_periods),dtype=float))
    p.linear_constraints.add(lin_expr=[
            cplex.SparsePair(flatten(B_soc_names[t]+
                    G_names[t]+P_minus_names_gen[t]+
                    P_PV_names[t]), val)
            for t in range(num_periods)], 
            senses=['G' for t in range(num_periods)], 
            rhs=scipy.zeros(num_periods), 
            names=['Sys_Op_15b_'+str(t+start_period)
                for t in range(num_periods)])
     
    #constraint 15c -- break purchasing symmetry
    for tech in technologies: 
        if tech.GetType() != "PVArray":
            for i in range(tech_maxes[tech.GetName()]-1): 
                ind = ["W_" + tech.GetName() + ".K" + str(i+1),
                    "W_" + tech.GetName() + ".K" + str(i+2)]
                val = [1,-1]
                p.linear_constraints.add(lin_expr=[
                    cplex.SparsePair( ind,val )],
                    senses = ['G'],
                    rhs = [0],
                    names = ['Sys_Op_15c_'+tech.GetName()+
                            ".K"+str(i)])
         

    ### Generator Operations               
    #constraint 16a -- Generator Max Power
    for gen in generators:
        for k in range(tech_maxes[gen.GetName()]):
            p.linear_constraints.add(lin_expr=[
                    cplex.SparsePair([G_var(gen,k+1,t+start_period),
                    P_minus_var(gen,k+1,(t+start_period))], 
                    [-1*gen.GetMaxPower(),1])
                    for t in range(num_periods)], 
            senses=['L' for t in range(num_periods)], 
            rhs=scipy.zeros(num_periods,dtype=float), 
            names=['Gen_Op_16a1_'+gen.GetName()+".K"
                +str(k+1)+".TT"+str(t+start_period)
                for t in range(num_periods)])
    
    #constraint 16a2 -- Generator Min Power
    for gen in generators:
        for k in range(tech_maxes[gen.GetName()]):
            p.linear_constraints.add(lin_expr=[
                    cplex.SparsePair([G_var(gen,k+1,t+start_period),
                    P_minus_var(gen,k+1,(t+start_period))], 
                    [-1*gen.GetMinPower(),1])
                    for t in range(num_periods)], 
            senses=['G' for t in range(num_periods)], 
            rhs=scipy.zeros(num_periods,dtype=float), 
            names=['Gen_Op_16a2_'+gen.GetName()+".K"
                +str(k+1)+".TT"+str(t+start_period)
                for t in range(num_periods)])
                    
    #constraint 16b -- Fuel Usage
    val = [1]
    val.extend(flatten([[gen.GetFuelUseB()/-1000.0 for i in range(tech_maxes[gen.GetName()])] 
            for gen in generators] ) )
    val.extend(flatten([[-1.0*gen.GetFuelUseC() for i in range(tech_maxes[gen.GetName()])] 
            for gen in generators] )) 
    
    p.linear_constraints.add(lin_expr=[
                    cplex.SparsePair([F_tilde_names[t]]+ 
                    flatten(P_minus_names_gen[t]
                    + G_names[t] ), val)
                    for t in range(num_periods)], 
            senses=['G' for t in range(num_periods)], 
            rhs=scipy.zeros(num_periods,dtype=float), 
            names=['Gen_Op_16b_'+str(t+start_period)
                for t in range(num_periods)])  
                
    #constraint 16c -- Use only generators that are purchased
    for gen in generators:
        for k in range(tech_maxes[gen.GetName()]):
            p.linear_constraints.add(lin_expr=[
                    cplex.SparsePair([
                    "G_" + gen.GetName() 
                    + ".K" + str(k+1) 
                    + ".TT" + str(t+start_period), "W_" + 
                    gen.GetName() + ".K" + str(k+1)],[1,-1])
                    for t in range(num_periods)],
                    senses=['L' for t in range(num_periods)], 
                rhs=scipy.zeros(num_periods,dtype=float), 
                names=['Gen_Op_16c_'+gen.GetName()+".K"
                +str(k+1)+".TT"+str(t+start_period)
                for t in range(num_periods)]) 
    
    #constraint 16d -- Break generator power out symmetry
    for gen in generators:
        for k in range(tech_maxes[gen.GetName()]-1):
            p.linear_constraints.add(lin_expr=[
                    cplex.SparsePair([
                    "P_minus_" + gen.GetName() + ".K" + str(k+1)
                    + ".TT" + str(t+start_period),
                    "P_minus_" + gen.GetName() 
                    + ".K" + str(k+2) 
                    + ".TT" + str(t+start_period)],[1,-1])
                    for t in range(num_periods)],
                    senses=['G' for t in range(num_periods)], 
                rhs=scipy.zeros(num_periods,dtype=float), 
                names=['Gen_Op_16d_'+gen.GetName()+".K"
                +str(k+2)+".TT"+str(t+start_period)
                for t in range(num_periods)])
    
    #constraint 16e -- Break generator use symmetry                                                                                                                                                      
    for gen in generators:
        for k in range(tech_maxes[gen.GetName()]-1):
            p.linear_constraints.add(lin_expr=[
                    cplex.SparsePair([
                    "G_" + gen.GetName() + ".K" + str(k+1)
                    + ".TT" + str(t+start_period),
                    "G_" + gen.GetName() 
                    + ".K" + str(k+2) 
                    + ".TT" + str(t+start_period)],[1,-1])
                    for t in range(num_periods)],
                    senses=['G' for t in range(num_periods)], 
                rhs=scipy.zeros(num_periods,dtype=float), 
                names=['Gen_Op_16d_'+gen.GetName()+".K"
                +str(k+2)+".TT"+str(t+start_period)
                for t in range(num_periods)])
    
    #constraint 17a -- PV out limited by solar power and # purchased
    for pv in pv_arrays:
        p.linear_constraints.add(lin_expr=[
                cplex.SparsePair(["P_PV_"+pv.GetName() 
                    + ".TT"+str(t+start_period),"X_"+pv.GetName()],
                    [1.0,
                    -1.0*pv.CalculatePowerAvailable(t+start_period-1)]) 
                    for t in range(num_periods)],
                senses=['L' for t in range(num_periods)],
                rhs=scipy.zeros(num_periods,dtype=float), 
                names=['PV_Op_17a_'+pv.GetName()
                    + ".TT"+str(t+start_period)
                    for t in range(num_periods)]
            )
    
    #constraint 17b -- upper bound on purchases
    for pv in pv_arrays:
        p.linear_constraints.add(lin_expr=[
                cplex.SparsePair(["X_"+pv.GetName()],
                    [1.0])],
                senses=['L'],
                rhs=[75],#[float(scalars["solar_max"])], 
                names=['PV_Op_17b_'+pv.GetName()]
            )
    
    #constraints 18a, 18b use Z in the place of B_soc*I_plus, substituted via
    #constraints 23-24, which are not implemented here for efficiency.
    #constraint Linear_18a: Linearization of P_plus (Z replaces B_soc*I_plus)
    for bat in batteries:
        val = [1, 
                -1.0*bat.GetVoltageA(),
                -1.0*bat.GetVoltageB()-(bat.GetAvgCurrent()*
                        bat.GetResistance())]
        for k in range(tech_maxes[bat.GetName()]):
            p.linear_constraints.add(
                    lin_expr=[cplex.SparsePair([
                            P_plus_var(bat,k+1,t),
                            Z_plus_var(bat,k+1,t),
                            I_plus_var(bat,k+1,t)]
                            ,val) for t in range(start_period,
                                end_period+1)],
                    senses=['E' for t in range(num_periods)], 
                    rhs=scipy.zeros(num_periods), 
                    names=['Linear_18a_'+bat.GetName()+".K"+str(k+1)
                        +".TT"+str(t) for t in range(start_period,
                                end_period+1)])
    #constraint Linear_18b: Linearization of P_minus (Z replaces B_soc*I_minus)
    for bat in batteries:
        val = [1, 
                -1.0*bat.GetVoltageA(),
                -1.0*bat.GetVoltageB()+(bat.GetAvgCurrent()*
                        bat.GetResistance())]
        for k in range(tech_maxes[bat.GetName()]):
            p.linear_constraints.add(
                    lin_expr=[cplex.SparsePair([
                            P_minus_var(bat,k+1,t),
                            Z_minus_var(bat,k+1,t),
                            I_minus_var(bat,k+1,t)]
                            ,val) for t in range(start_period,
                                end_period+1)],
                    senses=['E' for t in range(num_periods)], 
                    rhs=scipy.zeros(num_periods), 
                    names=['Linear_2_'+bat.GetName()+".K"+str(k+1)
                        +".TT"+str(t) for t in range(start_period,
                                end_period+1)])
    
    #constraint 18c -- State of charge tracking
    for bat in batteries:
        val = [1.0,-1.0,-1.0*(
            float(scalars["tau"])*bat.GetEfficiencyIn()
            / bat.GetReferenceCapacity()),
            1.0*(float(scalars["tau"]) / bat.GetReferenceCapacity())]
        for k in range(tech_maxes[bat.GetName()]):
            p.linear_constraints.add(lin_expr=[
                    cplex.SparsePair([B_soc_var(bat,k+1,t),
                                B_soc_var(bat,k+1,t-1),
                                I_plus_var(bat,k+1,t),
                                I_minus_var(bat,k+1,t)]
                            ,val) for t in range(start_period,end_period+1)],
                    senses=['E' for t in range(num_periods)],
                    rhs=scipy.zeros(num_periods), 
                    names=['Bat_Op_18c_'+bat.GetName()+".K"+str(k+1) +
                ".TT"+str(t) for t in range(start_period,end_period+1)])
                    
    for bat in batteries:
        #constraint 18d1 -- SOC maximum (note: includes soc_start)
        val = [1,-1.0*float(scalars["soc_max"])]
        for k in range(tech_maxes[bat.GetName()]):
            p.linear_constraints.add(lin_expr=[
                    cplex.SparsePair([B_soc_var(bat,k+1,t),
                                W_var(bat,k+1)]
                            ,val) for t in range(start_period-1,end_period+1)],
                    senses=['L' for t in range(num_periods+1)],
                    rhs=scipy.zeros(num_periods+1), 
                    names=['Bat_Op_18d1_'+bat.GetName()+".K"+str(k+1) +
                ".TT"+str(t) for t in range(start_period-1,end_period+1)])
    for bat in batteries:
        #constraint 18d2 -- SOC minimum (note: includes soc_start)
        val = [1,-1.0*float(scalars["soc_min"])]
        for k in range(tech_maxes[bat.GetName()]):
            p.linear_constraints.add(lin_expr=[
                    cplex.SparsePair([B_soc_var(bat,k+1,t),
                                W_var(bat,k+1)]
                            ,val) for t in range(start_period-1,end_period+1)],
                    senses=['G' for t in range(num_periods+1)],
                    rhs=scipy.zeros(num_periods+1), 
                    names=['Bat_Op_18d2_'+bat.GetName()+".K"+str(k+1) +
                ".TT"+str(t) for t in range(start_period-1,end_period+1)])
    
    for bat in batteries:
        #constraint 18e -- SOC equal for all batteries of type b
        val = [1.0,-1.0,1.0]
        for k in range(tech_maxes[bat.GetName()]-1):
            p.linear_constraints.add(lin_expr=[
                    cplex.SparsePair([B_soc_var(bat,k+1,t),
                                B_soc_var(bat,k+2,t),
                                W_var(bat,k+1)]
                        ,val) for t in range(start_period,end_period+1)],
                    senses=['L' for t in range(num_periods+1)],
                    rhs=1.0*scipy.ones(num_periods), 
                    names=['Bat_Op_18e_'+bat.GetName()+".K"+str(k+1) +
                ".TT"+str(t) for t in range(start_period,end_period+1)])
    for bat in batteries:
        #constraint 18f -- SOC equal for all batteries of type b
        val = [1.0,1.0,-1.0]
        for k in range(tech_maxes[bat.GetName()]-1):
            p.linear_constraints.add(lin_expr=[
                    cplex.SparsePair([B_soc_var(bat,k+1,t),
                                B_soc_var(bat,k+2,t),
                                W_var(bat,k+1)]
                        ,val) for t in range(start_period,end_period+1)],
                    senses=['G' for t in range(num_periods)],
                    rhs=-1.0*scipy.ones(num_periods), 
                    names=['Bat_Op_18f_'+bat.GetName()+".K"+str(k+1) +
                ".TT"+str(t) for t in range(start_period,end_period+1)])
 
    #constraint 18g -- tracks starting storage inventory
    ind = ['battery_inv']  #R^ell
    val = [-1.0]          
    for bat in batteries:
        for k in range(1,tech_maxes[bat.GetName()]+1):
            ind.append(soc_start_var(bat,k))
            val.append(bat.GetReferenceCapacity())
    p.linear_constraints.add(
        lin_expr=[cplex.SparsePair(ind,val)],
        senses=['E'], 
        rhs=[0], 
        names=['INVENTORY_18g'])
    
    #constraint 18h -- tracks ending storage inventory
    ind = ['battery_inv']
    val = [-1.0]          
    for bat in batteries:
        for k in range(1,tech_maxes[bat.GetName()]+1):
            ind.append(soc_end_var(bat,k))
            val.append(bat.GetReferenceCapacity())
    p.linear_constraints.add(
        lin_expr=[cplex.SparsePair(ind,val)],
        senses=['E'], 
        rhs=[0], 
        names=['INVENTORY_18h'])
    
    for bat in batteries:
        #constraint 18i1 -- power out upper limits of battery type b
        val = [1.0,-1.0*bat.GetMaxPower()]
        for k in range(tech_maxes[bat.GetName()]):
            p.linear_constraints.add(lin_expr=[
                    cplex.SparsePair([P_plus_var(bat,k+1,t),
                                B_plus_var(bat,k+1,t)]
                        ,val) for t in range(start_period,end_period+1)],
                    senses=['L' for t in range(num_periods)],
                    rhs=scipy.zeros(num_periods), 
                    names=['Bat_Op_18i1_'+bat.GetName()+".K"+str(k+1) +
                ".TT"+str(t) for t in range(start_period,end_period+1)])
    for bat in batteries:
        #constraint 18i2 -- power in lower limits of battery type b
        val = [1.0,-1.0*bat.GetMinPower()]
        for k in range(tech_maxes[bat.GetName()]):
            p.linear_constraints.add(lin_expr=[
                    cplex.SparsePair([P_plus_var(bat,k+1,t),
                                B_plus_var(bat,k+1,t)]
                        ,val) for t in range(start_period,end_period+1)],
                    senses=['G' for t in range(num_periods)],
                    rhs=scipy.zeros(num_periods), 
                    names=['Bat_Op_18i2_'+bat.GetName()+".K"+str(k+1) +
                ".TT"+str(t) for t in range(start_period,end_period+1)])
                    
    for bat in batteries:
        #constraint 18j1 -- power in upper limits of battery type b
        val = [1.0,-1.0*bat.GetMaxPower()]
        for k in range(tech_maxes[bat.GetName()]):
            p.linear_constraints.add(lin_expr=[
                    cplex.SparsePair([P_minus_var(bat,k+1,t),
                                B_minus_var(bat,k+1,t)]
                        ,val) for t in range(start_period,end_period+1)],
                    senses=['L' for t in range(num_periods)],
                    rhs=scipy.zeros(num_periods), 
                    names=['Bat_Op_18j1_'+bat.GetName()+".K"+str(k+1) +
                ".TT"+str(t) for t in range(start_period,end_period+1)])
    for bat in batteries:
        #constraint 18j2 -- power out lower limits of battery type b
        val = [1.0,-1.0*bat.GetMinPower()]
        for k in range(tech_maxes[bat.GetName()]):
            p.linear_constraints.add(lin_expr=[
                    cplex.SparsePair([P_minus_var(bat,k+1,t),
                                B_minus_var(bat,k+1,t)]
                        ,val) for t in range(start_period,end_period+1)],
                    senses=['G' for t in range(num_periods)],
                    rhs=scipy.zeros(num_periods), 
                    names=['Bat_Op_18j2_'+bat.GetName()+".K"+str(k+1) +
                ".TT"+str(t) for t in range(start_period,end_period+1)])
                    
                    
    for bat in batteries:
        #constraint 18k1 -- current in upper limits (semicontinuous)
        val = [1.0,-1.0*bat.GetIUPlus()]
        for k in range(tech_maxes[bat.GetName()]):
            p.linear_constraints.add(lin_expr=[
                    cplex.SparsePair([I_plus_var(bat,k+1,t),
                                B_plus_var(bat,k+1,t)]
                            ,val) for t in range(start_period,end_period+1)],
                    senses=['L' for t in range(num_periods)],
                    rhs=scipy.zeros(num_periods), 
                    names=['Bat_Op_18k1_'+bat.GetName()+".K"+str(k+1) +
                ".TT"+str(t) for t in range(start_period,end_period+1)])
    for bat in batteries:
        #constraint 18k2 -- current in lower limits (semicontinuous)
        val = [1.0,-1.0*bat.GetMinCurrent()]
        for k in range(tech_maxes[bat.GetName()]):
            p.linear_constraints.add(lin_expr=[
                    cplex.SparsePair([I_plus_var(bat,k+1,t),
                                B_plus_var(bat,k+1,t)]
                            ,val) for t in range(start_period,end_period+1)],
                    senses=['G' for t in range(num_periods)],
                    rhs=scipy.zeros(num_periods), 
                    names=['Bat_Op_18k2_'+bat.GetName()+".K"+str(k+1) +
                ".TT"+str(t) for t in range(start_period,end_period+1)])
                    
    for bat in batteries:
        #constraint 18l1 -- current out upper limits (semicontinuous)
        val = [1.0,-1.0*bat.GetIUMinus()]
        for k in range(tech_maxes[bat.GetName()]):
            p.linear_constraints.add(lin_expr=[
                    cplex.SparsePair([I_minus_var(bat,k+1,t),
                                B_minus_var(bat,k+1,t)]
                            ,val) for t in range(start_period,end_period+1)],
                    senses=['L' for t in range(num_periods)],
                    rhs=scipy.zeros(num_periods), 
                    names=['Bat_Op_18l1_'+bat.GetName()+".K"+str(k+1) +
                ".TT"+str(t) for t in range(start_period,end_period+1)])
    for bat in batteries:
        #constraint 18l2 -- current out upper limits (semicontinuous)
        val = [1.0,-1.0*bat.GetMinCurrent()]
        for k in range(tech_maxes[bat.GetName()]):
            p.linear_constraints.add(lin_expr=[
                    cplex.SparsePair([I_minus_var(bat,k+1,t),
                                B_minus_var(bat,k+1,t)]
                            ,val) for t in range(start_period,end_period+1)],
                    senses=['G' for t in range(num_periods)],
                    rhs=scipy.zeros(num_periods), 
                    names=['Bat_Op_18l2_'+bat.GetName()+".K"+str(k+1) +
                ".TT"+str(t) for t in range(start_period,end_period+1)])

    for bat in batteries:
        #constraint 18m -- current out upper limits (bound by SOC)
        val = [1.0,-1.0*(bat.GetIUMinus())]
        for k in range(tech_maxes[bat.GetName()]):
            p.linear_constraints.add(lin_expr=[
                    cplex.SparsePair([I_minus_var(bat,k+1,t),
                                B_soc_var(bat,k+1,t-1)]
                        ,val) for t in range(start_period,end_period+1)],
                    senses=['L' for t in range(num_periods)],
                    rhs=scipy.zeros(num_periods), 
                    names=['Bat_Op_18m_'+bat.GetName()+".K"+str(k+1) +
                ".TT"+str(t) for t in range(start_period,end_period+1)])
    
    for bat in batteries:
        #constraint 18n -- Only use batteries that you buy
        val = [1,1,-1]
        for k in range(tech_maxes[bat.GetName()]):
            p.linear_constraints.add(lin_expr=[
                    cplex.SparsePair([B_minus_var(bat,k+1,t),
                                B_plus_var(bat,k+1,t),
                                W_var(bat,k+1)],val) 
                        for t in range(start_period,end_period+1)],
                    senses=['L' for t in range(num_periods)],
                    rhs=scipy.zeros(num_periods), 
                    names=['Bat_Op_18n_'+bat.GetName()+".K"+str(k+1) +
                ".TT"+str(t) for t in range(start_period,end_period+1)])
        
    for bat in batteries:
        #constraint 18o -- Can't simultaneously charge and discharge batteries 
        val=[1,1]
        for batp in batteries:
            if batp != bat:
                for k in range(tech_maxes[bat.GetName()]):
                    for kp in range(tech_maxes[batp.GetName()]):
                        p.linear_constraints.add(lin_expr=[
                            cplex.SparsePair([B_plus_var(bat,k+1,t),
                            B_minus_var(batp,kp+1,t)],val) 
                            for t in range(
                                start_period,end_period+1)],
                    senses=['L' for t in range(num_periods)],
                    rhs=scipy.ones(num_periods,dtype=float), 
                    names=['Bat_Op_18o_'+bat.GetName()+".K"+str(k+1) 
                        +".BP"+batp.GetName()[-1]+
                        ".KP"+str(kp+1)+".TT"+str(t) 
                        for t in range(start_period,end_period+1)])
    
    
    #Constraint 19a: Nonanticipativity of W      
    #Fix purchase decision if running the upper bound model
    if upper_bound:
        inds = []
        vals = []
        for tech in technologies: 
            if tech.GetType() != "PVArray":
                for k in range(1,tech_maxes[tech.GetName()]+1):
                    inds.append(W_var(tech,k))
                    if tech.GetName() in W_init[scenario].keys(): 
                        if "K"+str(k) in W_init[scenario][tech.GetName()].keys():
                            vals.append(W_init[scenario][tech.GetName()]["K"+str(k)])
                        else: vals.append(0)
                    else: vals.append(0)
        #print inds, vals
        p.variables.set_upper_bounds(zip(inds,vals))
        p.variables.set_lower_bounds(zip(inds,vals))  
    
    #Constraint 19b: Nonanticipativity of X   
    #Fix purchase decision if running the upper bound model
    if upper_bound:
        inds = []
        vals = []
        for tech in technologies: 
            if tech.GetType() == "PVArray":
                inds.append(X_var(tech))
                if tech.GetName() in X_init[scenario].keys():
                    vals.append(X_init[scenario][tech.GetName()])
                else: vals.append(0) 
        #print inds, vals
        p.variables.set_upper_bounds(zip(inds,vals))
        p.variables.set_lower_bounds(zip(inds,vals))  
    
    #Constraint 19c: Nonanticipativity of R^\ell.  This is enforced by 
    #fixing the starting state of charge for all batteriesin the upper bound.
    #model, and is relaxed in the lower bound model.
    if upper_bound and start_period != 1:
        val = [1.0,-1.0*boundary_soc]
        sense = ['E']
        for bat in batteries:
            for k in range(1,tech_maxes[bat.GetName()]+1):
                p.linear_constraints.add(
                    lin_expr=[cplex.SparsePair([soc_start_var(bat,k),W_var(bat,k)],
                        val)],
                    senses=sense, 
                    rhs=[0], 
                    names=['Storage_inv_start_UB_19c1'+bat.GetName()+".K"+str(k)])
    if upper_bound and start_period != 1:
        val = [1.0,-1.0*boundary_soc]
        sense = ['E']
        ending = 'UB_'
        for bat in batteries:
            for k in range(1,tech_maxes[bat.GetName()]+1):
                p.linear_constraints.add(
                    lin_expr=[cplex.SparsePair([soc_end_var(bat,k),W_var(bat,k)],
                        val)],
                    senses=sense, 
                    rhs=[0], 
                    names=['Storage_inv_end_UB_19c2'+bat.GetName()+".K"+str(k)]) 
   
    #Constraint 20: boundary condition for battery state-of-charge
    if start_period == 1: 
        val = [1.0,-1.0*scalars["b_init"]]
        sense = ['E']
        for bat in batteries:
            for k in range(1,tech_maxes[bat.GetName()]+1):
                p.linear_constraints.add(
                    lin_expr=[cplex.SparsePair([soc_start_var(bat,k),W_var(bat,k)],
                        val)],
                    senses=sense, 
                    rhs=[0], 
                    names=['Storage_inv_boundary_'+bat.GetName()+".K"+str(k)])  
        
    #constraints (21) denote variable bounds
    #constraints (22)-(24) not directly implemented    
    
    def i_u_part_minus(bat,n):
        """
        This defines:
        i^{L-}_{bn} = i_u_part_minus(bat,n-1), and
        i^{U-}_{bn} = i_u_part_minus(n)
        """
        return bat.GetMinCurrent() + batPartsMinus[n]*(bat.GetIUMinus()-bat.GetMinCurrent())
    def i_u_part_plus(bat,n):
        """
        This defines:
        i^{L+}_{bn} = i_u_part_plus(bat,n-1), and
        i^{U+}_{bn} = i_u_part_plus(n)
        """
        return bat.GetMinCurrent() + batPartsPlus[n]*(bat.GetIUPlus()-bat.GetMinCurrent())
    
    
    #constraint 25a -- one lambda per BKT - reconcile plus.
    val = [1.0 for i in range(1,numBatParts+1)]+[-1.0]
    for bat in batteries:
        for k in range(1,tech_maxes[bat.GetName()]+1):
             p.linear_constraints.add(
                    lin_expr=[cplex.SparsePair([lambda_plus_var(bat,k,n,t) 
                        for n in range(1,numBatParts+1)]+[B_plus_var(bat,k,t)],
                            val) for t in range(start_period,
                                    end_period+1)],
                    senses=['E' for t in range(num_periods)], 
                    rhs=scipy.zeros(num_periods), 
                    names=['Linear_25a_'+bat.GetName()+".K"+str(k)+
                        ".TT"+str(t) for t in range(start_period,
                                    end_period+1)])       
            
    #constraint 25b -- one lambda per BKT - reconcile minus.
    for bat in batteries:
        for k in range(1,tech_maxes[bat.GetName()]+1):
             p.linear_constraints.add(
                    lin_expr=[cplex.SparsePair([lambda_minus_var(bat,k,n,t) 
                        for n in range(1,numBatParts+1)]+[B_minus_var(bat,k,t)],
                            val) for t in range(start_period,
                                    end_period+1)],
                    senses=['E' for t in range(num_periods)], 
                    rhs=scipy.zeros(num_periods), 
                    names=['Linear_25b_'+bat.GetName()+".K"+str(k)+
                        ".TT"+str(t) for t in range(start_period,
                                    end_period+1)]) 
    
    #constraint 25c -- lambda_plus assignment based on I_plus (lower bounds).
    for bat in batteries:
        for n in range(1,numBatParts+1):
            val = [1.0,
                -1.0*bat.GetMinCurrent(),
                -1.0*(i_u_part_plus(bat,n-1)-bat.GetMinCurrent())]
            for k in range(1,tech_maxes[bat.GetName()]+1):
                p.linear_constraints.add(
                    lin_expr=[cplex.SparsePair([
                            I_plus_var(bat,k,t),
                            B_plus_var(bat,k,t),
                            lambda_plus_var(bat,k,n,t)],
                            val) for t in range(start_period,
                                    end_period+1)],
                    senses=['G' for t in range(num_periods)], 
                    rhs=scipy.zeros(num_periods), 
                    names=['Linear_25c_'+bat.GetName()+".K"+str(k)+
                        ".N"+str(n)+
                        ".TT"+str(t) for t in range(start_period,
                                    end_period+1)])
                        
    #constraint 25d -- lambda_plus assignment based on I_plus (upper bounds).
    for bat in batteries:
        for n in range(1,numBatParts+1):
            val = [1.0,
                -1.0*bat.GetIUPlus(),
                1.0*(bat.GetIUPlus()-i_u_part_plus(bat,n))]
            for k in range(1,tech_maxes[bat.GetName()]+1):
                p.linear_constraints.add(
                    lin_expr=[cplex.SparsePair([
                            I_plus_var(bat,k,t),
                            B_plus_var(bat,k,t),
                            lambda_plus_var(bat,k,n,t)],
                            val) for t in range(start_period,
                                    end_period+1)],
                    senses=['L' for t in range(num_periods)], 
                    rhs=scipy.zeros(num_periods), 
                    names=['Linear_25d_'+bat.GetName()+".K"+str(k)+
                        ".N"+str(n)+
                        ".TT"+str(t) for t in range(start_period,
                                    end_period+1)])
                        
    #constraint 25e -- lambda_minus assignment based on I_minus (lower bounds)
    for bat in batteries:
        for n in range(1,numBatParts+1):
            val = [1.0,
                -1.0*bat.GetMinCurrent(),
                -1.0*(i_u_part_minus(bat,n-1)-bat.GetMinCurrent())]
            for k in range(1,tech_maxes[bat.GetName()]+1):
                p.linear_constraints.add(
                    lin_expr=[cplex.SparsePair([
                            I_minus_var(bat,k,t),
                            B_minus_var(bat,k,t),
                            lambda_minus_var(bat,k,n,t)],
                            val) for t in range(start_period,
                                    end_period+1)],
                    senses=['G' for t in range(num_periods)], 
                    rhs=scipy.zeros(num_periods), 
                    names=['Linear_25e_'+bat.GetName()+".K"+str(k)+
                        ".N"+str(n)+
                        ".TT"+str(t) for t in range(start_period,
                                    end_period+1)])
    
    #constraint 25f -- lambda_minus assignment based on I_minus (upper bounds)
    for bat in batteries:
        for n in range(1,numBatParts+1):
            val = [1.0,
                -1.0*bat.GetIUMinus(),
                1.0*(bat.GetIUMinus()-i_u_part_minus(bat,n))]
            for k in range(1,tech_maxes[bat.GetName()]+1):
                p.linear_constraints.add(
                    lin_expr=[cplex.SparsePair([
                            I_minus_var(bat,k,t),
                            B_minus_var(bat,k,t),
                            lambda_minus_var(bat,k,n,t)],
                            val) for t in range(start_period,
                                    end_period+1)],
                    senses=['L' for t in range(num_periods)], 
                    rhs=scipy.zeros(num_periods), 
                    names=['Linear_25f_'+bat.GetName()+".K"+str(k)+
                        ".N"+str(n)+
                        ".TT"+str(t) for t in range(start_period,
                                    end_period+1)])
                                   
    ##constraint 25g -- Z_plus linearization
    for bat in batteries:
        for n in range(1,numBatParts+1):
            val = [1.0,
                -1.0*i_u_part_plus(bat,n),
                1.0*bat.GetMaxSOC(),
                -1.0*(bat.GetMaxSOC()-bat.GetMinSOC())*
                    (bat.GetIUPlus()-bat.GetMinCurrent()),
                1.0*(bat.GetMaxSOC()-bat.GetMinSOC())*
                    (bat.GetIUPlus()-bat.GetMinCurrent()) + 
                bat.GetMaxSOC()*i_u_part_plus(bat,n)
                ]
            for k in range(1,tech_maxes[bat.GetName()]+1):
                p.linear_constraints.add(
                    lin_expr=[cplex.SparsePair([
                            Z_plus_var(bat,k,t),
                            B_soc_var(bat,k,t-1),
                            I_plus_var(bat,k,t),
                            lambda_plus_var(bat,k,n,t),
                            W_var(bat,k)],
                            val) for t in range(start_period,
                                    end_period+1)],
                    senses=['G' for t in range(num_periods)], 
                    rhs=scipy.zeros(num_periods), 
                    names=['Linear_25g_'+bat.GetName()+".K"+str(k)+
                        ".N"+str(n)+
                        ".TT"+str(t) for t in range(start_period,
                                    end_period+1)])
    ##constraint 25h-- Z_plus linearization
    for bat in batteries:
        for n in range(1,numBatParts+1):
            val = [1.0,
                -1.0*i_u_part_plus(bat,n-1),
                -1.0*bat.GetMinSOC(),
                -1.0*(bat.GetMaxSOC()-bat.GetMinSOC())*
                    (bat.GetIUPlus()-bat.GetMinCurrent()),
                1.0*(bat.GetMaxSOC()-bat.GetMinSOC())*
                    (bat.GetIUPlus()-bat.GetMinCurrent())+
                    bat.GetMinSOC()*i_u_part_plus(bat,n-1)
                ]
            for k in range(1,tech_maxes[bat.GetName()]+1):
                p.linear_constraints.add(
                    lin_expr=[cplex.SparsePair([
                            Z_plus_var(bat,k,t),
                            B_soc_var(bat,k,t-1),
                            I_plus_var(bat,k,t),
                            lambda_plus_var(bat,k,n,t),
                            W_var(bat,k)],
                            val) for t in range(start_period,
                                    end_period+1)],
                    senses=['G' for t in range(num_periods)], 
                    rhs=scipy.zeros(num_periods), 
                    names=['Linear_25h_'+bat.GetName()+".K"+str(k)+
                        ".N"+str(n)+
                        ".TT"+str(t) for t in range(start_period,
                                    end_period+1)])  
    ##constraint 25i -- Z_plus linearization
    for bat in batteries:
        for n in range(1,numBatParts+1):
            val = [1.0,
                -1.0*i_u_part_plus(bat,n-1),
                -1.0*bat.GetMaxSOC(),
                1.0*(bat.GetMaxSOC()-bat.GetMinSOC())*
                    (bat.GetIUPlus()-bat.GetMinCurrent()),
                -1.0*(bat.GetMaxSOC()-bat.GetMinSOC())*
                    (bat.GetIUPlus()-bat.GetMinCurrent())+
                    bat.GetMaxSOC()*i_u_part_plus(bat,n-1)
                ]
            for k in range(1,tech_maxes[bat.GetName()]+1):
                p.linear_constraints.add(
                    lin_expr=[cplex.SparsePair([
                            Z_plus_var(bat,k,t),
                            B_soc_var(bat,k,t-1),
                            I_plus_var(bat,k,t),
                            lambda_plus_var(bat,k,n,t),
                            W_var(bat,k)],
                            val) for t in range(start_period,
                                    end_period+1)],
                    senses=['L' for t in range(num_periods)], 
                    rhs=scipy.zeros(num_periods), 
                    names=['Linear_25i_'+bat.GetName()+".K"+str(k)+
                        ".N"+str(n)+
                        ".TT"+str(t) for t in range(start_period,
                                    end_period+1)])
    ##constraint 25j -- Z_plus linearization
    for bat in batteries:
        for n in range(1,numBatParts+1):
            val = [1.0,
                -1.0*i_u_part_plus(bat,n),
                -1.0*bat.GetMinSOC(),
                1.0*(bat.GetMaxSOC()-bat.GetMinSOC())*
                    (bat.GetIUPlus()-bat.GetMinCurrent()),
                -1.0*(bat.GetMaxSOC()-bat.GetMinSOC())*
                    (bat.GetIUPlus()-bat.GetMinCurrent())+
                    bat.GetMinSOC()*i_u_part_plus(bat,n)
                ]
            for k in range(1,tech_maxes[bat.GetName()]+1):
                p.linear_constraints.add(
                    lin_expr=[cplex.SparsePair([
                            Z_plus_var(bat,k,t),
                            B_soc_var(bat,k,t-1),
                            I_plus_var(bat,k,t),
                            lambda_plus_var(bat,k,n,t),
                            W_var(bat,k)],
                            val) for t in range(start_period,
                                    end_period+1)],
                    senses=['L' for t in range(num_periods)], 
                    rhs=scipy.zeros(num_periods), 
                    names=['Linear_25j_'+bat.GetName()+".K"+str(k)+
                        ".N"+str(n)+
                        ".TT"+str(t) for t in range(start_period,
                                    end_period+1)])   
    ##constraint 25k -- Z_minus linearization
    for bat in batteries:
        for n in range(1,numBatParts+1):
            val = [1.0,
                -1.0*i_u_part_minus(bat,n),
                1.0*bat.GetMaxSOC(),
                -1.0*(bat.GetMaxSOC()-bat.GetMinSOC())*
                    (bat.GetIUMinus()-bat.GetMinCurrent()),
                1.0*(bat.GetMaxSOC()-bat.GetMinSOC())*
                    (bat.GetIUMinus()-bat.GetMinCurrent())+
                    bat.GetMaxSOC()*i_u_part_minus(bat,n)
                ]
            for k in range(1,tech_maxes[bat.GetName()]+1):
                p.linear_constraints.add(
                    lin_expr=[cplex.SparsePair([
                            Z_minus_var(bat,k,t),
                            B_soc_var(bat,k,t-1),
                            I_minus_var(bat,k,t),
                            lambda_minus_var(bat,k,n,t),
                            W_var(bat,k)],
                            val) for t in range(start_period,
                                    end_period+1)],
                    senses=['G' for t in range(num_periods)], 
                    rhs=scipy.zeros(num_periods), 
                    names=['Linear_25k_'+bat.GetName()+".K"+str(k)+
                        ".N"+str(n)+
                        ".TT"+str(t) for t in range(start_period,
                                    end_period+1)])
    ##constraint 25l-- Z_minus linearization
    for bat in batteries:
        for n in range(1,numBatParts+1):
            val = [1.0,
                -1.0*i_u_part_minus(bat,n-1),
                -1.0*bat.GetMinSOC(),
                -1.0*(bat.GetMaxSOC()-bat.GetMinSOC())*
                    (bat.GetIUMinus()-bat.GetMinCurrent()),
                1.0*(bat.GetMaxSOC()-bat.GetMinSOC())*
                    (bat.GetIUMinus()-bat.GetMinCurrent())+
                            bat.GetMinSOC()*i_u_part_minus(bat,n-1)        
                ]
            for k in range(1,tech_maxes[bat.GetName()]+1):
                p.linear_constraints.add(
                    lin_expr=[cplex.SparsePair([
                            Z_minus_var(bat,k,t),
                            B_soc_var(bat,k,t-1),
                            I_minus_var(bat,k,t),
                            lambda_minus_var(bat,k,n,t),
                            W_var(bat,k)],
                            val) for t in range(start_period,
                                    end_period+1)],
                    senses=['G' for t in range(num_periods)], 
                    rhs=scipy.zeros(num_periods),
                    names=['Linear_25l_'+bat.GetName()+".K"+str(k)+
                        ".N"+str(n)+
                        ".TT"+str(t) for t in range(start_period,
                                    end_period+1)])  
    ##constraint 25m -- Z_minus linearization
    for bat in batteries:
        for n in range(1,numBatParts+1):
            val = [1.0,
                -1.0*i_u_part_minus(bat,n-1),
                -1.0*bat.GetMaxSOC(),
                1.0*(bat.GetMaxSOC()-bat.GetMinSOC())*
                    (bat.GetIUMinus()-bat.GetMinCurrent()),
                -1.0*(bat.GetMaxSOC()-bat.GetMinSOC())*
                    (bat.GetIUMinus()-bat.GetMinCurrent())+
                    bat.GetMaxSOC()*i_u_part_minus(bat,n-1)        
                ]
            for k in range(1,tech_maxes[bat.GetName()]+1):
                p.linear_constraints.add(
                    lin_expr=[cplex.SparsePair([
                            Z_minus_var(bat,k,t),
                            B_soc_var(bat,k,t-1),
                            I_minus_var(bat,k,t),
                            lambda_minus_var(bat,k,n,t),
                            W_var(bat,k)],
                            val) for t in range(start_period,
                                    end_period+1)],
                    senses=['L' for t in range(num_periods)], 
                    rhs=scipy.zeros(num_periods), 
                    names=['Linear_25m_'+bat.GetName()+".K"+str(k)+
                        ".N"+str(n)+
                        ".TT"+str(t) for t in range(start_period,
                                    end_period+1)])
    ##constraint 25n-- Z_minus linearization
    for bat in batteries:
        for n in range(1,numBatParts+1):
            val = [1.0,
                -1.0*i_u_part_minus(bat,n),
                -1.0*bat.GetMinSOC(),
                1.0*(bat.GetMaxSOC()-bat.GetMinSOC())*
                    (bat.GetIUMinus()-bat.GetMinCurrent()),
                -1.0*(bat.GetMaxSOC()-bat.GetMinSOC())*
                    (bat.GetIUMinus()-bat.GetMinCurrent())+
                    bat.GetMinSOC()*i_u_part_minus(bat,n)  
                ]
            for k in range(1,tech_maxes[bat.GetName()]+1):
                p.linear_constraints.add(
                    lin_expr=[cplex.SparsePair([
                            Z_minus_var(bat,k,t),
                            B_soc_var(bat,k,t-1),
                            I_minus_var(bat,k,t),
                            lambda_minus_var(bat,k,n,t),
                            W_var(bat,k)],
                            val) for t in range(start_period,
                                    end_period+1)],
                    senses=['L' for t in range(num_periods)], 
                    rhs=scipy.zeros(num_periods), 
                    names=['Linear_25n_'+bat.GetName()+".K"+str(k)+
                        ".N"+str(n)+
                        ".TT"+str(t) for t in range(start_period,
                                    end_period+1)])   
    
    #constraint 25o -- Z_minus <= I_minus
    for bat in batteries:
        val = [1.0,-1.0]  
        for k in range(1,tech_maxes[bat.GetName()]+1):
            p.linear_constraints.add(
                lin_expr=[cplex.SparsePair([
                        Z_minus_var(bat,k,t),
                        I_minus_var(bat,k,t)],
                        val) for t in range(start_period,
                                end_period+1)],
                senses=['L' for t in range(num_periods)], 
                rhs=scipy.zeros(num_periods), 
                names=['Linear_25o_'+bat.GetName()+".K"+str(k)+
                    ".TT"+str(t) for t in range(start_period,
                                end_period+1)])    
    
    #constraint 25p -- Z_plus <= I_plus    
    for bat in batteries:
        val = [1.0,-1.0]  
        for k in range(1,tech_maxes[bat.GetName()]+1):
            p.linear_constraints.add(
                lin_expr=[cplex.SparsePair([
                        Z_plus_var(bat,k,t),
                        I_plus_var(bat,k,t)],
                        val) for t in range(start_period,
                                end_period+1)],
                senses=['L' for t in range(num_periods)], 
                rhs=scipy.zeros(num_periods), 
                names=['Linear_25p_'+bat.GetName()+".K"+str(k)+
                    ".TT"+str(t) for t in range(start_period,
                                end_period+1)])             
                
    
    #constraint CUT1 -- lower bounds on number of generators purchased;
    #see Scioletti et al., 2017, Section 3.1, for how this cut is generated
    ind = []
    val = []
    for gen in generators:
        for k in range(1,tech_maxes[gen.GetName()]+1):
            ind.append(W_var(gen,k)) 
            val.append(1) 
    p.linear_constraints.add(
                lin_expr=[cplex.SparsePair(ind,val)],
                senses=['G'], 
                rhs=[mingens], 
                #rhs = [0],
                names=['CUT1'])
    #constraint CUT2 -- at most 1 battery is allowed (taken from Scioletti et al. 2017, OPTE)
    ind = []
    val = []
    for bat in batteries:
        for k in range(1,tech_maxes[bat.GetName()]+1):
            ind.append(W_var(bat,k)) 
            val.append(1) 
    p.linear_constraints.add(
                lin_expr=[cplex.SparsePair(ind,val)],
                senses=['L'], 
                rhs=[1], 
                names=['CUT2'])   
    #Constraint RESET_SOCa - fix SOC at beginning of period or restrict to purchased batteries
    if start_period == 1: 
        val = [1.0,-1.0*scalars["b_init"]]
        sense = ['E']
        ending = 'START_'
    elif upper_bound:
        val = [1.0,-1.0*boundary_soc]
        sense = ['E']
        ending = 'UB_'
    else:
        val = [1.0,-1.0]
        sense = ['L']
        ending = 'LB_'
    for bat in batteries:
        for k in range(1,tech_maxes[bat.GetName()]+1):
            p.linear_constraints.add(
                lin_expr=[cplex.SparsePair([soc_start_var(bat,k),W_var(bat,k)],
                    val)],
                senses=sense, 
                rhs=[0], 
                names=['Storage_inv_start_'+ending+bat.GetName()+".K"+str(k)])  
    #Constraint RESET_SOCb - fix SOC at end of period
    if upper_bound:
        val = [1.0,-1.0*boundary_soc]
        sense = ['E']
        ending = 'UB_'
    else:
        val = [1.0,-1.0]
        sense = ['L']
        ending = 'LB_' 
    for bat in batteries:
        for k in range(1,tech_maxes[bat.GetName()]+1):
            p.linear_constraints.add(
                lin_expr=[cplex.SparsePair([soc_end_var(bat,k),W_var(bat,k)],
                    val)],
                senses=sense, 
                rhs=[0], 
                names=['RESET_SOC_end_'+ending+bat.GetName()+".K"+str(k)]) 
                                
    #Valid Inequality: Constraint RESET_POINTS_EQ - start_soc and end_soc must be equal
    if start_period != 1:
        val = [1.0,-1.0]
        sense = ['E']
        for bat in batteries:
            for k in range(1,tech_maxes[bat.GetName()]+1):
                p.linear_constraints.add(
                    lin_expr=[cplex.SparsePair([soc_end_var(bat,k),soc_start_var(bat,k)],
                        val)],
                    senses=sense, 
                    rhs=[0], 
                    names=['START_END_SOC_EQ_'+bat.GetName()+".K"+str(k)])
    
    #constraint INVENTORY -- tracks ending storage inventory (which is the same
    #as starting inventory)
    ind = ['battery_inv']
    val = [-1.0]          
    for bat in batteries:
        for k in range(1,tech_maxes[bat.GetName()]+1):
            ind.append(soc_end_var(bat,k))
            val.append(bat.GetReferenceCapacity())
    p.linear_constraints.add(
        lin_expr=[cplex.SparsePair(ind,val)],
        senses=['E'], 
        rhs=[0], 
        names=['INVENTORY'])
    
    #Constraint GEN_CAPACITY - set minimum generator capacity; See Section
    #4.4.1 of the paper for how this cut is generated.  Implementation is 
    #noted in the module 'RunPH'
    ind = []
    val = []
    for gen in generators:
        for k in range(1,tech_maxes[gen.GetName()]+1):
            ind.append(W_var(gen,k))
            val.append(gen.GetMaxPower())
    p.linear_constraints.add(
        lin_expr = [cplex.SparsePair(ind,val)],
        senses = ['G'],
        #rhs = [0],
        rhs = [mincap],
        names = ['GEN_MIN_CAP']
        )
        
    #Fix purchase decision if running the upper bound formulation
    if upper_bound:
        inds = []
        vals = []
        for tech in technologies: 
            if tech.GetType() != "PVArray":
                for k in range(1,tech_maxes[tech.GetName()]+1):
                    inds.append(W_var(tech,k))
                    if tech.GetName() in W_init[scenario].keys(): 
                        if "K"+str(k) in W_init[scenario][tech.GetName()].keys():
                            vals.append(W_init[scenario][tech.GetName()]["K"+str(k)])
                        else: vals.append(0)
                    else: vals.append(0)
            else: 
                inds.append(X_var(tech))
                if tech.GetName() in X_init[scenario].keys():
                    vals.append(X_init[scenario][tech.GetName()])
                else: vals.append(0) 
        #print inds, vals
        p.variables.set_upper_bounds(zip(inds,vals))
        p.variables.set_lower_bounds(zip(inds,vals))  
    
    #If there is a candidate design decision given as input, we test the 
    #feasibility by fixing the solution and observing the CPLEX output.
    inds = []
    vals = []
    if IFS != None:
        for k in IFS.keys():
            if type(IFS[k]) == list:
                for ind,val in IFS[k]:
                    inds.append(ind)
                    vals.append(val)
        p.MIP_starts.add(cplex.SparsePair(inds,vals), p.MIP_starts.effort_level.solve_MIP)
        
        
    # SETTING BRANCHING PRIORITIES
    #p.order.set([(W_names[i],1,p.order.branch_direction.up) for i in range(len(W_names))])
    
    # Specify the problem type and solve
    p.set_problem_type( p.problem_type.MILP )
    if output: p.write("stochmodelcur.lp")
    sol = {}
    clock = time.time() 
    try: 
        p.solve()
        elapsed = time.time() - clock
        #output all variables in solution  
        sol["start"]=start_period
        sol["end"]=end_period
        sol["W"] = [(name,int(0.5+p.solution.get_values(name))) for name in flatten(W_names)]
        sol["obj"] = p.solution.get_objective_value() 
        sol["lb"] = p.solution.MIP.get_best_objective()
        sol["gap"] = p.solution.MIP.get_mip_relative_gap()
        sol["X"] = [(name,p.solution.get_values(name)) for name in X_names]
        sol["P_minus"] = [(name,p.solution.get_values(name)) for name in flatten(P_minus_names) if p.solution.get_values(name) > 1e-6]
        sol["P_plus"] = [(name,p.solution.get_values(name)) for name in flatten(P_plus_names) if p.solution.get_values(name) > 1e-6]
        sol["P_PV"] = [(name,p.solution.get_values(name)) for name in flatten(P_PV_names) if p.solution.get_values(name) > 1e-6]
        sol["F_tilde"] = [(name,p.solution.get_values(name)) for name in F_tilde_names if p.solution.get_values(name) > 1e-7]
        sol["B_soc"] = [(name,p.solution.get_values(name)) for name in flatten(B_soc_names) if p.solution.get_values(name) > 1e-8]
        sol["I_minus"] = [(name,p.solution.get_values(name)) for name in flatten(I_minus_names) if p.solution.get_values(name) > 1e-7]
        sol["I_plus"] = [(name,p.solution.get_values(name)) for name in flatten(I_plus_names) if p.solution.get_values(name) > 1e-7]
        sol["Z_minus"] = [(name,p.solution.get_values(name)) for name in flatten(Z_minus_names) if p.solution.get_values(name) > 1e-7]
        sol["Z_plus"] = [(name,p.solution.get_values(name)) for name in flatten(Z_plus_names) if p.solution.get_values(name) > 1e-7]
        sol["B_minus"] = [(name,p.solution.get_values(name)) for name in flatten(B_minus_names) if p.solution.get_values(name) > 0.5]
        sol["B_plus"] = [(name,p.solution.get_values(name)) for name in flatten(B_plus_names) if p.solution.get_values(name) > 0.5]
        sol["G"] = [(name,p.solution.get_values(name)) for name in flatten(G_names) if p.solution.get_values(name) > 0.5]
        sol["lambda_plus"] = [(name,p.solution.get_values(name)) for name in flatten(lambda_plus_names) if p.solution.get_values(name) > 0.5]
        sol["lambda_minus"] = [(name,p.solution.get_values(name)) for name in flatten(lambda_minus_names) if p.solution.get_values(name) > 0.5]
        sol["B_soc_start"] =sum([p.solution.get_values(name) for name in flatten(soc_start_names)])
        sol["B_soc_end"] = sum([p.solution.get_values(name) for name in flatten(soc_end_names)])
        sol["battery_inv"] = p.solution.get_values('battery_inv')
        sol["lamb_start"] = inv_mult
        sol["lamb_end"] = inv_mult
        sol["I_plus_max"] = sum([bat.GetIUPlus() for bat in batteries if W_var(bat,1) in flatten(W_names) and p.solution.get_values(W_var(bat,1)) > 0.5])
        sol["I_minus_max"] = sum([bat.GetIUMinus() for bat in batteries if W_var(bat,1) in flatten(W_names) and p.solution.get_values(W_var(bat,1)) > 0.5])
        sol["IL"] = sum([bat.GetMinCurrent() for bat in batteries if W_var(bat,1) in flatten(W_names) and p.solution.get_values(W_var(bat,1)) > 0.5])
    except cplex.exceptions.CplexSolverError:
        elapsed = time.time() - clock
        sol["W"] = None
        sol["obj"] = scipy.infty
        sol["lb"] = scipy.infty
    sol["solve_time"] = elapsed
    del p
    return sol     
    

def OutputUBSolutions(ub_sols,output_dir="GAMS/",scenario_name = "ll12"):
    W_file = open(output_dir+'WW.csv','w')
    X_file = open(output_dir+'X.csv', 'w')
#    L_file = open(output_dir+'L.csv', 'w')
    sols = sorted([(sol["start"],sol) for sol in ub_sols])
    append = False
#    L_vals = {}
    X_vals = sols[0][1]["X"]
    W_vals = sols[0][1]["W"]
    for sol in sols:
#        x = sol[1]
#        for name, val in x:
#            if name in L_vals.keys():
#                L_vals[name] += val
#            else:
#                L_vals[name] = val
        OutputSolution(sol,append,output_dir,scenario_name)
        append = True
    import re
    def RemoveVarNameHeader(var_name,scenario_name = "LL12"):
        """Returns the variable indices without the variable name header."""
        return re.findall(r"^.*\_([A-Z0-9\.]*)",var_name)[-1]
    for name, val in W_vals:
        W_file.write(scenario_name+"."+RemoveVarNameHeader(name)+"  "+str(val)+"\n")
    for name, val in X_vals:
        X_file.write(scenario_name+"."+RemoveVarNameHeader(name)+"  "+str(val)+"\n")
#    for name in L_vals.keys():
#        L_file.write(scenario_name+"."+RemoveVarNameHeader(name)+"  "+str(L_vals[name])+"\n")
        
    
def OutputSolution(sol, append=False, output_dir = "./../sol_files/", scenario_name = "ll12"):
    """Writes variable values to files.
    sol -- solution taken from CPLEX solver
    append -- indicator of whether to append to file or write to a new one
    output_lead -- directory link to output files destination
    scenario_name -- name of scenario, for testing in GAMS
    retval -- None, output written to files"""
    if append: app = 'a'
    else: app = 'w'
    P_minus_file = open(output_dir+'P_out.csv',app) 
    P_solar_file = open(output_dir+'P_solar.csv',app) 
    P_plus_file = open(output_dir+'P_in.csv',app)
    F_tilde_file = open(output_dir+'F_tilde.csv',app)
    B_soc_file = open(output_dir+'B_soc.csv',app)
    B_plus_file = open(output_dir+'B_plus.csv',app)
    B_minus_file = open(output_dir+'B_minus.csv',app)
    G_file = open(output_dir+'G.csv',app)
    I_plus_file = open(output_dir+'I_plus.csv',app)
    I_minus_file = open(output_dir+'I_minus.csv',app)
    Z_plus_file = open(output_dir+'Z_plus.csv',app)
    Z_minus_file = open(output_dir+'Z_minus.csv',app)
#    Y_plus_file = open(output_dir+'Y_plus.csv',app)
#    Y_minus_file = open(output_dir+'Y_minus.csv',app)
    #lamb_plus_file = open(output_dir+'lamb_plus.csv',app)
    #lamb_minus_file = open(output_dir+'lamb_minus.csv',app)
    #Z_file = open(output_dir+'Z.csv',app)
    #Z_file.write(scenario_name + "  " + str(sol["obj"])+"\n")
    
    import re
    def RemoveVarNameHeader(var_name,scenario_name = "LL12"):
        """Returns the variable indices without the variable name header."""
        return re.findall(r"^.*\_([A-Z0-9\.]*)",var_name)[-1]
    for name, val in sol["P_minus"]:
        P_minus_file.write(scenario_name+"."+RemoveVarNameHeader(name)+"  "+str(val)+"\n")
    for name, val in sol["P_plus"]:
        P_plus_file.write(scenario_name+"."+RemoveVarNameHeader(name)+"  "+str(val)+"\n")
    for name, val in sol["P_PV"]:
        P_solar_file.write(scenario_name+"."+RemoveVarNameHeader(name)+"  "+str(val)+"\n")
    for name, val in sol["F_tilde"]:
        F_tilde_file.write(scenario_name+"."+RemoveVarNameHeader(name)+"  "+str(val)+"\n")
    for name, val in sol["B_soc"]:
        B_soc_file.write(scenario_name+"."+RemoveVarNameHeader(name)+"  "+str(val)+"\n")
    for name, val in sol["B_minus"]:
        B_minus_file.write(scenario_name+"."+RemoveVarNameHeader(name)+"  "+str(val)+"\n")
    for name, val in sol["B_plus"]:
        B_plus_file.write(scenario_name+"."+RemoveVarNameHeader(name)+"  "+str(val)+"\n")
    for name, val in sol["I_minus"]:
        I_minus_file.write(scenario_name+"."+RemoveVarNameHeader(name)+"  "+str(val)+"\n")
    for name, val in sol["I_plus"]:
        I_plus_file.write(scenario_name+"."+RemoveVarNameHeader(name)+"  "+str(val)+"\n")
    #for name, val in sol["V_soc"]:
    #    V_soc_file.write(scenario_name+"."+RemoveVarNameHeader(name)+"  "+str(val)+"\n")
#    for name, val in sol["I_minus"]:
#        Y_minus_file.write(scenario_name+"."+RemoveVarNameHeader(name)+"  "+str(val)+"\n")
#    for name, val in sol["I_plus"]:
#        Y_plus_file.write(scenario_name+"."+RemoveVarNameHeader(name)+"  "+str(val)+"\n")
    for name, val in sol["Z_minus"]:
        Z_minus_file.write(scenario_name+"."+RemoveVarNameHeader(name)+"  "+str(val)+"\n")
    for name, val in sol["Z_plus"]:
        Z_plus_file.write(scenario_name+"."+RemoveVarNameHeader(name)+"  "+str(val)+"\n")
    for name, val in sol["G"]:
        G_file.write(scenario_name+"."+RemoveVarNameHeader(name)+"  "+str(val)+"\n")
    #for name, val in sol["lambda_plus"]:
    #    lamb_plus_file.write(scenario_name+"."+RemoveVarNameHeader(name)+"  "+str(val)+"\n")
    #for name, val in sol["lambda_minus"]:
    #    lamb_minus_file.write(scenario_name+"."+RemoveVarNameHeader(name)+"  "+str(val)+"\n")
        
def OutputAllSols(sols, output_dir = "GAMS/", scenario_name = "ll12"):
    """Writes variable values to files, for all solutions.
    sols -- list of solutions taken from CPLEX solver
    output_lead -- directory link to output files destination
    scenario_name -- name of scenario, for testing in GAMS
    retval -- None, output written to files"""
    Z = 0.0
#    L = {}
    W = {}
    import re
    def RemoveVarNameHeader(var_name,scenario_name = "LL12"):
        """Returns the variable indices without the variable name header."""
        return re.findall(r"^.*\_([A-Z0-9\.]*)",var_name)[-1]
    for name, val in sols[0]["W"]: 
        #print name
#        L[RemoveVarNameHeader(name)] = 0.0
        W[RemoveVarNameHeader(name)] = val
    append = False
    for sol in sols:
        #print sol["start"]
        OutputSolution(sol,append=append, output_dir=output_dir, scenario_name=scenario_name)
        Z += sol["obj"]
#        for name, val in sol["L"]:
#            L[RemoveVarNameHeader(name)] += val
        append = True
    Z_file = open(output_dir+'Z.csv','w')
    Z_file.write(scenario_name + "  " + str(Z)+"\n")
#    L_file = open(output_dir+'L.csv','w')
    W_file = open(output_dir+'WW.csv','w')
    X_file = open(output_dir+'X.csv','w')
    for name in W.keys():
#        L_file.write(scenario_name+"."+name+"  "+str(L[name])+"\n")
        W_file.write(scenario_name+"."+name+"  "+str(W[name])+"\n")
    for name, val in sols[0]["X"]:
        X_file.write(scenario_name+"."+RemoveVarNameHeader(name)+"  "+str(val)+"\n")
    Z_file.close()
    W_file.close()
    X_file.close()

    
def SolveScenarioLBPartition(params):  
    scalar_filename = "./../inputs/Opt_Model_Scalars.csv"
    load_filename = "./../inputs/Opt_Model_Demand.csv"
    fc_filename = "./../inputs/Opt_Model_Fuel.csv"
    tech_filename = "./../inputs/Opt_Model_J_inputs.csv"
    pv_filename = "./../inputs/Opt_Model_Solar.csv"
    pv_spec_filename = "./../inputs/Opt_Model_PV_Inputs.csv"
    max_filename = "./../inputs/Opt_Model_J_MAX.csv"  
    cuts_filename = "./../inputs/Opt_Model_Cuts.csv" 
    mincap_filename = "./../inputs/mingencap.csv" 
    inputs = InputReader.GetInputs(scalar_filename,load_filename,
            pv_filename,pv_spec_filename,fc_filename,tech_filename,
            max_filename,cuts_filename,8760) 
    InputReader.setBounds(inputs,params[4])
    sol = SolvePartition(inputs,params[4],
        params[0],min(params[1],8760),params[2],batPartsPlus=params[5],
        batPartsMinus=params[6],design_mults=params[7],upper_bound=False, 
        mipgap=0.005,reslim=params[10],IFS=params[8],mincap=params[9]
        )
    #print "LB", sol["lb"]
    return sol
    
def SolveScenarioUBPartition(params): 
    if local:  
        #filenames for running locally
        scalar_filename = "./../inputs/Opt_Model_Scalars.csv"
        load_filename = "./../inputs/Opt_Model_Demand.csv"
        fc_filename = "./../inputs/Opt_Model_Fuel.csv"
        tech_filename = "./../inputs/Opt_Model_J_inputs.csv"
        pv_filename = "./../inputs/Opt_Model_Solar.csv"
        pv_spec_filename = "./../inputs/Opt_Model_PV_Inputs.csv"
        max_filename = "./../inputs/Opt_Model_J_MAX.csv"  
        cuts_filename = "./../inputs/Opt_Model_Cuts.csv" 
        W_filename = "./../sol_files/"+params[4]+"/WW.csv" 
        X_filename = "./../sol_files/"+params[4]+"/X.csv" 
    else:
        #filenames for running on HPC resources
        scalar_filename = "./../inputs/Opt_Model_Scalars.csv"
        load_filename = "./../inputs/Opt_Model_Demand.csv"
        fc_filename = "./../inputs/Opt_Model_Fuel.csv"
        tech_filename = "./../inputs/Opt_Model_J_inputs.csv"
        pv_filename = "./../inputs/Opt_Model_Solar.csv"
        pv_spec_filename = "./../inputs/Opt_Model_PV_Inputs.csv"
        max_filename = "./../inputs/Opt_Model_J_MAX.csv"  
        cuts_filename = "./../inputs/Opt_Model_Cuts.csv" 
        W_filename = "./../inputs/sol_files/"+params[4]+"/WW.csv" 
        X_filename = "./../inputs/sol_files/"+params[4]+"/X.csv" 
        mincap_filename = "./../inputs/mingencap.csv" 
    inputs = InputReader.GetInputs(scalar_filename,load_filename,
            pv_filename,pv_spec_filename,fc_filename,tech_filename,
            max_filename,cuts_filename,8760, W_filename=W_filename,
            X_filename=X_filename)
    #print X_filename,inputs["W_init"],inputs["X_init"] 
    InputReader.setBounds(inputs,params[4])
    sol = SolvePartition(inputs,params[4],
        params[0],min(params[1],8760),params[2],batPartsPlus=params[5],
        batPartsMinus=params[6],upper_bound=True,mipgap=0.005,
        reslim=params[8], boundary_soc = params[7])
    return sol
         
   
def RunLinearPartitionsWeek(scenario):
    timeclock = time.time()
    #set running parameters
    params_list = [(i*168+1,(i+1)*168,0.0,scenario) for i in range(3)]
    params_list.append((8737,8760,0.0,scenario))
    obj=0.0
    #append=False
    for param_tuple in params_list:
        print (param_tuple)
        sol = SolveScenarioUBPartition(param_tuple)
        #OutputSolution(sol,append)
        obj += sol["obj"]
        #append=True
    print ("Final UB Objective Value:",obj) 
    ub =obj
    ubt = time.time() - timeclock 
    timeclock = time.time()
    obj=0.0
    print ("Lower Bound calculating:")
    for param_tuple in params_list:
        print (param_tuple)
        sol = SolveScenarioLBPartition(param_tuple)
        #OutputSolution(sol,append)
        obj += sol["obj"]
        #append=True
    lb = obj
    lbt = time.time() - timeclock
    print ("Final LB Objective Value:",obj)
    print ("LB Time:", lbt)
    print ("UB Time:", ubt)
    print ("Total Time: " + str(ubt+lbt))
    print ("Gap: "+str((ub-lb)/ub)) 
    
def RunParallelLBDay(scenario, lambs=[0.0]*365, 
            procs=multiprocessing.cpu_count(),num_days=365, probSize=24,
             batPartsPlus=[0,0.25,0.5,0.75,1.0],batPartsMinus=[0,0.25,0.5,0.75,1.0],
             socParts=[0.0,0.5,1.0],
             design_mults=[(0,0)]*365,amortize_on_demand = True, IFS=None, mincap=0,
             reslim=60):
    if IFS == None: IFS = [None]*num_days
    params_list = [(i*probSize+1,(i+1)*probSize,lambs[i],0.0,scenario,
            batPartsPlus,batPartsMinus,design_mults[i],IFS[i],mincap,reslim,socParts) for i in range(num_days)]
    pool = multiprocessing.Pool(processes=procs)
    results = pool.map(SolveScenarioLBPartition, params_list)
    lb = sum(sol['lb'] for sol in results)
    frac_over = sum([1.0 for sol in results if sol["solve_time"] >= reslim]) / num_days
    pool.close()
    pool.terminate()
    del pool
    sortedResults = sorted([(x["start"],x) for x in results])
    results = [x[1] for x in sortedResults]
    return lb, results, frac_over
    
def RunParallelUBDay(scenario, lambs=[0.0]*365,  procs=multiprocessing.cpu_count(), 
        incumbent=scipy.infty,num_days=365, probSize=24, 
        batPartsPlus=[0,0.25,0.5,0.75,1.0],batPartsMinus=[0,0.25,0.5,0.75,1.0],
        socParts=[0.0,0.5,1.0],
        boundary_soc=0.5,reslim=60):
    params_list = [(i*probSize+1,(i+1)*probSize,lambs[i],0.0,scenario,
            batPartsPlus,batPartsMinus,boundary_soc,reslim,socParts) for i in range(num_days)]
    pool = multiprocessing.Pool(processes=procs)
    #t = time.time()
    results = pool.map(SolveScenarioUBPartition, params_list)
    obj = sum(sol['obj'] for sol in results)
    frac_over = sum([1.0 for sol in results if sol["solve_time"] >= reslim]) / num_days
    pool.close()
    pool.terminate()
    del pool
    sortedResults = sorted([(x["start"],x) for x in results])
    results = [x[1] for x in sortedResults]
    return obj, results, frac_over  

def RunLinearLBWeek(scenario,lamb=0.0,batPartsPlus=[0,0.25,0.5,0.75,1.0],batPartsMinus=[0,0.25,0.5,0.75,1.0]):
    #params_list = [(i*168+1,(i+1)*168,lamb) for i in range(4)]
    params_list = [(i*168+1,(i+1)*168,lamb,scenario,batPartsPlus,batPartsMinus) for i in range(3)]
    params_list.append((8737,8760,lamb,scenario))
    obj = 0.0
    lb = 0.0
    for param_tuple in params_list:
        #print param_tuple
        sol = SolveScenarioLBPartition(param_tuple)
        #OutputSolution(sol,append)
        obj += sol["obj"]
        lb += sol["lb"]
        #append=True
    OutputSolution(sol)
    print ("lambda",params_list[0][2],"obj",obj,"lb",lb)
    return lb, sol
    
def RunLinearLBDay(scenario, lambs=[0.0]*365, 
        procs=multiprocessing.cpu_count(),num_days=365, probSize=24,
        batPartsPlus=[0,0.25,0.5,0.75,1.0],batPartsMinus=[0,0.25,0.5,0.75,1.0],
        design_mults=[(0,0)]*365,amortize_on_demand = True, IFS=None, mincap=0,
        reslim=180):
    if IFS == None: IFS = [None]*num_days
    params_list = [(i*probSize+1,(i+1)*probSize,lambs[i],0.0,scenario,
            batPartsPlus,batPartsMinus,design_mults[i],IFS[i],mincap,reslim) for i in range(num_days)]
    obj = 0.0
    lb = 0.0
    for param_tuple in params_list:
        #print param_tuple
        sol = SolveScenarioLBPartition(param_tuple)
        #OutputSolution(sol,append)
        obj += sol["obj"]
        lb += sol["lb"]
        print(sol["W"])
        #append=True
    if lb < scipy.infty: 
        OutputSolution(sol,scenario_name=scenario)
    print ("lambda",params_list[0][2],"obj",obj,"lb",lb)
    return lb, sol
    
def RunLinearUBDay(scenario,num_days=365,probSize=24,lamb=0.0,incumbent = scipy.infty,batPartsPlus=[0,0.25,0.5,0.75,1.0],batPartsMinus=[0,0.25,0.5,0.75,1.0]):
    params_list = [(i*probSize+1,(i+1)*probSize,lamb,lamb,scenario,batPartsPlus,batPartsMinus) for i in range(num_days)]
    #params_list = [(i*24+1,(i+1)*24,lamb,scenario) for i in [1,2,363,364]]
    obj = 0.0
    lb = 0.0
    for param_tuple in params_list:
        #print param_tuple
        sol = SolveScenarioUBPartition(param_tuple)
        #OutputSolution(sol,append)
        obj += sol["obj"]
        lb += sol["lb"]
        print(sol["W"])
        #append=True
    if lb < scipy.infty: 
        OutputSolution(sol,scenario_name=scenario)
    print ("lambda",params_list[0][2],"obj",obj,"lb",lb)
    return obj, sol
    
def RunLBBisection(scenario,lambrange=100.0,lambtol=1):
    import time 
    t = time.time()
    dist = lambrange/4.0
    curlamb = lambrange/2.0
    obj_lo,sol_lo = RunParallelLBDay(scenario,lamb=0.0)
    obj_hi,sol_hi = RunParallelLBDay(scenario,lamb=lambrange)
    obj_mid, sol_mid = RunParallelLBDay(scenario,curlamb)
    #obj_lo,sol_lo = RunLinearLBDay(scenario,lamb=0.0)
    #obj_hi,sol_hi = RunLinearLBDay(scenario,lamb=lambrange)
    #obj_mid, sol_mid = RunLinearLBDay(scenario,curlamb)
    if obj_mid > obj_lo and obj_mid > obj_hi:
        incumbent = sol_mid
        obj = obj_mid
        best_lam = curlamb
    elif obj_hi > obj_lo:
        incumbent = sol_hi
        obj = obj_hi
        best_lam = lambrange
    else: 
        incumbent = sol_lo
        obj = obj_lo
        best_lam = 0.0
    while dist >= lambtol:
        if obj_lo < obj_hi and obj_mid < obj_hi:  #increase lambda
            obj = obj_hi
            incumbent = sol_hi
            best_lam = curlamb + 2*dist
            obj_lo = obj_mid
            curlamb += dist
            obj_mid, sol_mid = RunParallelLBDay(scenario,curlamb)
        elif obj_lo > obj_mid and obj_lo > obj_hi:  #decrease lambda
            obj = obj_lo
            incumbent = sol_lo
            best_lam = curlamb - 2*dist
            obj_hi = obj_mid
            curlamb -= dist
            obj_mid, sol_mid = RunParallelLBDay(scenario,curlamb)
        else:  #mid is best, look on either side
            obj = obj_mid
            incumbent = sol_mid
            best_lam = curlamb
            obj_lo, sol_lo = RunParallelLBDay(scenario,curlamb-dist)
            obj_hi, sol_hi = RunParallelLBDay(scenario,curlamb+dist)
        dist /= 2.0
    print ("Best obj:", obj)
    elapsed = time.time() - t
    print ("Best lambda:", best_lam)
    print ("Time elapsed:", elapsed)
    return obj, incumbent, best_lam
    
def OutputDesign(sol,output_dir,scenario_name):
    import re
    def RemoveVarNameHeader(var_name,scenario_name = "LL12"):
        """Returns the variable indices without the variable name header."""
        return re.findall(r"^.*\_([A-Z0-9\.]*)",var_name)[-1]
    W_file = open(output_dir+'WW.csv','w')
    X_file = open(output_dir+'X.csv','w')
    W = {}
    X = {}
    for name, val in sol[0]: 
        W[RemoveVarNameHeader(name)] = val
    for name, val in sol[1]:
        X[RemoveVarNameHeader(name)] = int(0.5+val)
    for name in W.keys():
        W_file.write(scenario_name+"."+name+"  "+str(W[name])+"\n")
    for name in X.keys():
        X_file.write(scenario_name+"."+name+"  "+str(X[name])+"\n")
    W_file.close()
    X_file.close()
     
def RunLBSubgradient(scenario,num_days=365,probSize=24,itermax=40,tol=0.05,n=3,
        ubCheck=5,batPartsPlus=[0,0.25,0.5,0.75,1.0],batPartsMinus=[0,0.25,0.5,0.75,1.0]):
    lambs = scipy.array([0.0]*(num_days+1))
    obj,sol_ub = RunParallelUBDay(scenario,lambs=lambs,num_days=num_days,
            probSize=probSize,batPartsPlus=batPartsPlus,batPartsMinus=batPartsMinus)
    outputSolsLambs(scenario,sol_ub,0,True)
    sol_ub_sorted = sorted([(x["start"],x) for x in sol_ub])
    ubs = scipy.array([x[1]["obj"] for x in sol_ub_sorted])
    #print ubs
    #print scipy.argmax(ubs)
    ub_sols = [str(sol_ub[0]["W"])+str(sol_ub[0]["X"])]
    print ('upper bound', obj)
    lb,sol_lb = RunParallelLBDay(scenario,lambs=lambs,num_days=num_days,
            batPartsPlus=batPartsPlus,batPartsMinus=batPartsMinus,IFS=sol_ub)
    outputSolsLambs(scenario,sol_lb,0,False)
    incumbent=lb*1.0
    gap = (obj-incumbent)/obj
    k=1  #iteration number
    it=0  # number of iterations since last improvement of bound
    t = 0.5 # step size
    while gap > tol and k <= itermax:
        #print "updating lambda"
        lambs = UpdateLambda(lambs,sol_lb,ubs,t)
        print ('iter',k)
        #print 'first 12 lambs:',lambs[:12]
        #del lb
        #del sol_lb 
        lb,sol_lb = RunParallelLBDay(scenario,lambs=lambs,num_days=num_days,
                probSize=probSize,batPartsPlus=batPartsPlus,batPartsMinus=batPartsMinus,IFS=sol_ub)
        outputSolsLambs(scenario,sol_lb,k,False)
        print ("lb",lb)
        if lb > incumbent:
            incumbent = lb*1.0
            it=0
        else:
            it += 1
        if it >= n:
            print ('halving timestep - no improvement.')
            t /= 2
            it = 0
        if k % ubCheck == 0:
            sol, freq = GetMostPopularDesign(sol_lb)
            if sol not in ub_sols:
                print ("attempting new candidate.")
                objU,solU = RunParallelUBDay(scenario,num_days=num_days,
                        probSize=probSize,batPartsPlus=batPartsPlus,batPartsMinus=batPartsMinus)
                ub_sols.append(str(solU[0]["W"])+str(solU[0]["X"]))
                if objU <  obj:
                    print ("new incumbent upper bound found.")
                    obj = objU
                    sol_ub = solU
        k += 1
        gap = (obj-incumbent)/obj
        print ('updated gap:', gap)
    #OutputUBSolutions(ub_sols)
    return incumbent,obj,gap,(k-1)

def outputSolsLambs(scenario,sols,k,ub):
    """output key attributed of solutions (lambdas, soc start/end, ub and lb) to files."""
    if ub: writetype = 'w'
    else: writetype = 'a'
    lbfile = open("Design"+scenario+"lb.csv",writetype)
    lambfile = open("Design"+scenario+"lambs.csv",writetype)
    stfile = open("Design"+scenario+"starts.csv",writetype)
    endfile = open("Design"+scenario+"ends.csv",writetype)
    gapfile = open("Design"+scenario+"gaps.csv",writetype)
    if ub:
        lambfile.write('iter')
        stfile.write('iter')
        endfile.write('iter')
        lbfile.write('iter')
        gapfile.write('iter')
        for i in range(len(sols)):
            lambfile.write(","+str(i+1))
            stfile.write(","+str(i+1))
            endfile.write(","+str(i+1))
            lbfile.write(","+str(i+1))
            gapfile.write(","+str(i+1))
        lambfile.write('\nUB,')
        stfile.write('\nUB,')
        endfile.write('\nUB,')
        lbfile.write('\nUB,')
        gapfile.write('\nUB,')
    else:
        lambfile.write(str(k)+",")
        stfile.write(str(k)+",")
        endfile.write(str(k)+",")
        lbfile.write(str(k)+",") 
        gapfile.write(str(k)+",") 
    for sol in sols:
        if ub: lambfile.write(str(sol["lamb_end"])+",")
        stfile.write(str(sol["B_soc_start"])+",")
        endfile.write(str(sol["B_soc_end"])+",")
        if ub: lbfile.write(str(sol["obj"])+",")
        else: lbfile.write(str(sol["lb"])+",")
        gapfile.write(str(sol["gap"])+",")
    lambfile.write('\n')
    stfile.write('\n')
    endfile.write('\n')
    lbfile.write('\n')
    gapfile.write('\n')
    lambfile.close()
    stfile.close()
    endfile.close()
    lbfile.close() 
    gapfile.close()       
                        
def GetMostPopularDesign(sols):
    """ Returns the most popular design decision when given a collection of
    solutions to the lower bound subproblems as input.
    sols -- solutions to lower bound subproblems
    retval -- a tuple of the binary/integer solutions and the top frequency"""
    solutions = {}
    for sol in sols:
        if str(sol["W"])+str(sol["X"]) in solutions.keys():
            solutions[str(sol["W"])+str(sol["X"])] += 1
        else:
            solutions[str(sol["W"])+str(sol["X"])] = 1
    incumbent = None
    hiFreq = 0
    for sol in solutions.keys():
        if solutions[sol] > hiFreq:
            hiFreq = solutions[sol]
            incumbent = sol
    return incumbent, hiFreq
        
def UpdateLambda(lambs,sol_lb,ubs,t):
    sol = sorted((x["start"],x) for x in sol_lb)
    soc_starts = scipy.array([x[1]["B_soc_start"] for x in sol[1:]])
    soc_ends = scipy.array([x[1]["B_soc_end"] for x in sol[:-1]])
    lbs = scipy.array([x[1]["lb"] for x in sol])
    stepsizes = t *(ubs[1:]+ubs[:-1]-lbs[1:]-lbs[:-1])*(soc_starts-soc_ends)
    new_lambs = scipy.maximum(scipy.zeros_like(soc_starts),scipy.array(lambs[1:-1]) + stepsizes)
    return [0.0] + list(new_lambs) + [0]
    
def RunWholeProblem(scenario):
    timeclock = time.time()
    params_list = [(1,8760,0.0,0.0,scenario,[0.0,1.0])]
    obj=0.0
    append=False
    for param_tuple in params_list:
        print (param_tuple)
        sol = SolveScenarioLBPartition(param_tuple)
        OutputSolution(sol,append)
        obj += sol["obj"]
        append=True
    print ("Final UB Objective Value:",obj )
    ub =obj
    ubt = time.time() - timeclock 
    timeclock = time.time()
    print ("Total Time:",ubt)
    print ("OBJ", ub)
       
                                                                             
if __name__ == "__main__":
    import time
    import sys
    if len(sys.argv) > 1: scenario = sys.argv[1]
    else: scenario = "ll1"
#    clock = time.time()
    #print "Scenario:",scenario
    
    #QUICK TEST
    #scen_numbers = range(1,15)
    #scenario_names = ["ll"+str(x) for x in scen_numbers]
    #f = open("revised_method.csv","w")
    #f.write("scen,lb,time\n")
    #for scenario in scenario_names:
    #    print scenario
    #    clock = time.time()
    #    lb,obj = RunLinearLBDay(scenario,num_days = 3)
    #    elapsed = time.time()-clock
    #    f.write(scenario+","+str(lb)+","+str(elapsed)+"\n")
    #    print "SCEN:",scenario,"LB:",lb,"TIME:",elapsed
    #f.close()
    
    #Solve whole problem
    lb,obj = RunLinearLBDay(scenario, lambs=[0.0],num_days=1, probSize=24,
        batPartsPlus=[0,0.25,0.5,0.75,1.0],batPartsMinus=[0,0.25,0.5,0.75,1.0],
        design_mults=[(0,0)]*1,amortize_on_demand = True, IFS=None, mincap=0,
        reslim=60)
    print (scenario)
    
    #print "SUBGRADIENT ATTEMPT"
    #lb,obj,gap,iters = RunLBSubgradient(scenario,num_days=365,probSize=24,itermax=30,tol=0.05,n=3,ubCheck=5)
    
    #num_days = 1
    #probSize = 24
    #print "Grabbing all upper bounds."
    #scen_numbers = range(12,13)
    #scenario_names = ["ll"+str(x) for x in scen_numbers]
    ##output_leads = ["./../inputs/sol_files/ll"+str(x)+"/" for x in scen_numbers]
    #output_leads = ["sol/"]
    #for i in range(len(scen_numbers)):
    #    obj,sols = RunParallelUBDay(scenario_names[i],num_days=num_days,probSize=probSize)
    #    print "scenario",scenario_names[i]," UB completed."
    #    obj,sol_lb = RunParallelLBDay(scenario_names[i],num_days=num_days,probSize=probSize,IFS=sols)
    #    #print "obj:",obj
    #    OutputAllSols(sols, output_dir = output_leads[i], scenario_name = scenario_names[i])
    #    print "LB",scenario_names[i],"completed."
    #elapsed = time.time() - clock
    #print "TOTAL TIME:", elapsed
    
    ######PRINTING OUTPUT######
    #outfile = open("CurPartModel_R"+scenario[2:]+".csv",'w')
    #outfile.write("Scenario,LB,UB,gap,iters,time\n")
    #outfile.write(scenario+","+str(lb)+","+str(obj)+","+str(gap)+","+str(iters)+","+str(elapsed)+","+"\n")
    #
    #outfile.close()
    #print "Solving Whole Problem:"
    #clock = time.time()
    #RunWholeProblem("ll12")
    #elapsed = time.time() - clock
    #print "Completed in "+str(elapsed)+" second
    
    