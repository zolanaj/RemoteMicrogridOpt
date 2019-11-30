"""
This module will import methods that create CPLEX instances of differnet types,
and run them in sequence in a way such that the solution of one model informs 
the problem statement of the second.  Specifically, we run the upper bound
instance of the linear model (a.k.a. V1 of the optimization model) to obtain 
information on dispatch, and then use that in the updated model with a partition
scheme that puts bounds on the auxiliary bilinear terms in an effort to tighten 
the relaxation.

"""
import Partition_Univariate
import Partition_Univariate_G
import Partition_Bivariate
import Partition_Nagarajan
import Min_capacity_model
import scipy
import MIPtoMINLP
import time


def GetCurrentParts(sols,direction="I_minus",numParts=5):
    """  Obtains nonuniform current partition boundaries, using a CPLEX solution
    as input.
    sols - CPLEX solutions to subproblems
    direction - "I_plus" for current in; "I_minus" for current out
    numParts - number of partitions
    
    retval - partition boundaries to be used for partition on current
    """
    currents = []
    ub_cur = direction+"_max"
    for sol in sols:
        currents.extend([(x[1]-sol["IL"])/(sol[ub_cur]-sol["IL"]) for x in sol[direction]])
    sortedCurrents = sorted(currents)
    #print direction,sortedCurrents
    if len(currents) == 0: return [0.0,0.25,0.5,0.75,1.0]
    indices = [len(currents)*i/(numParts-1) for i in range(numParts-1)]+[len(currents)-1]
    parts = scipy.array([sortedCurrents[i] for i in indices])
    parts[0] = 0.0
    parts[-1] = 1.0
#    print (direction, parts)
    return parts 
    
def GetMinCap(scenario,lambs,design_mults,num_days=365,probSize=24,bat_parts=[0.0,1.0]):
    lb,sol_lb = Min_capacity_model.RunParallelLBDay(scenario,
            lambs=lambs,
            design_mults=design_mults,num_days=num_days,probSize=probSize,
            batPartsPlus=bat_parts,batPartsMinus=bat_parts)
    mincap = max([sol["lb"] for sol in sol_lb])
    #print mincaps
    return sol_lb, mincap
        
def RunPH(scenario,num_days=365,probSize=24,itermax=40,tol=0.05,n=3,ubCheck=5,
        tech_filename = "./../inputs/Opt_Model_J_inputs.csv",
        max_filename = "./../inputs/Opt_Model_J_MAX.csv",
        pv_spec_filename = "./../inputs/Opt_Model_PV_Inputs.csv",
        lbFind = 5,numBatParts=5, numSocParts=5, MINLP=True, nonuniform_part=False, findgen=True,
        timeLimit = 3600, solsToCheck=10, method = "univariate"):
    """Runs the subgradient algorithm for our model.
    scenario - scenario number/location
    num_days - number of subproblems to solve 
    probSize - size of single subproblem in time periods (e.g., hours)
    itermax - stopping criterion for algorithm: number of iterations of 
            progressive hedging
    tol - soptting criterion: optimality gap
    n - number of iterations before halving scalar step size multiplier
    ubCheck - number of iterations at which we attempt to find a new upper 
            bound
    tech_filename - technology spec sheet location
    pv_spec_filename - pv systems specs file location
    lbFind - # of iterations before finding a lower bound (removing proximal 
            term)
    numBatParts -- number of current subregions for current variable
    numSocParts -- number of current subregions for battery state-of-charge 
            variable
    MINLP -- determines whether or not an MINLP gap is sought
    nonuniform_part -- calculates a nonuniform partitioning scheme if True
    findgen -- Establishes a minimim generator capacity if True (used for 
                    generating a constraint in day-long subproblems)
    timeLimit -- time-based stopping criterion for total algorithm
    solsToCheck -- number of UB runs to find solutions
    method -- partititioning model selected for approximating bilinear 
            relationships
    """
    iterfile = open("iter-method-"+method+"NOVI-scen"+scenario+"-m"+str(numSocParts)+"-n"+str(numBatParts)+".csv",'a')
    iterfile.write("iteration,LB,time\n")
    iterfile.close()
    IFS_dir = "./../sol_files/"+scenario+"/"
    inputs_dir = "./../inputs/"
    clock = time.time()
    lambs = [0.0]*(num_days)
    design_mults, rho, rho_bat = GenerateDesignMultipliers(scenario=scenario,
                tech_filename = tech_filename,
                max_filename = max_filename,
                pv_spec_filename = pv_spec_filename,
                num_subproblems = num_days,
                probSize=probSize,
                proximal=False)
    #obtain minimum generator capacity.
    bat_parts_plus = scipy.arange(numBatParts+1)*1.0/numBatParts
    bat_parts_minus = bat_parts_plus * 1.0
    socParts = scipy.arange(numSocParts+1)*1.0/numSocParts
    if findgen:
        sols, mincap = GetMinCap(scenario,lambs,design_mults,num_days=num_days,probSize=probSize,bat_parts=[0.0,1.0])
        print ("mincap:",mincap)
        if nonuniform_part:
            bat_parts_plus = GetCurrentParts(sols,"I_plus",numBatParts) 
            bat_parts_minus = GetCurrentParts(sols,"I_minus",numBatParts) 
    else:
        mincap = 0
    #obtain upper bound.
    obj = scipy.infty
    ub_sols = []
    sol_ub = None

    incumbent=1e-6
    gap = (obj-incumbent)/incumbent
    k=0 #iteration number
    frac_overs = []
    #it=0  # number of iterations since last improvement of bound
    #t = 0.5 # step size
    
    elapsed = time.time()-clock
    while gap > tol and k <= itermax and elapsed < timeLimit:
        if method == "nagarajan":  #Model (N)
            lb,sol_lb,frac_over = Partition_Nagarajan.RunParallelLBDay(scenario,
                lambs=lambs,design_mults=design_mults,num_days=num_days,
                probSize=probSize,batPartsPlus=bat_parts_plus,
                batPartsMinus=bat_parts_minus,mincap=mincap,socParts=socParts)
        elif method == "univariate":  #Model (U)
            lb,sol_lb,frac_over = Partition_Univariate.RunParallelLBDay(scenario,
                lambs=lambs,design_mults=design_mults,num_days=num_days,
                probSize=probSize,batPartsPlus=bat_parts_plus,
                batPartsMinus=bat_parts_minus,mincap=mincap)
        elif method == "gounaris":    #Model (G)
            lb,sol_lb,frac_over = Partition_Univariate_G.RunParallelLBDay(scenario,
                lambs=lambs,design_mults=design_mults,num_days=num_days,
                probSize=probSize,batPartsPlus=bat_parts_plus,
                batPartsMinus=bat_parts_minus,mincap=mincap)
        elif method == "bivariate":    #Model (U), bivariate implementation
            lb,sol_lb,frac_over = Partition_Bivariate.RunParallelLBDay(scenario,
                lambs=lambs,design_mults=design_mults,num_days=num_days,
                probSize=probSize,batPartsPlus=bat_parts_plus,
                batPartsMinus=bat_parts_minus,mincap=mincap,socParts=socParts)
        else: 
            print("invalid method!")
            return -1
        elapsed = time.time() - clock
        frac_overs.append(frac_over)
        if nonuniform_part:
            bat_parts_plus = GetCurrentParts(sol_lb,"I_plus",numBatParts) 
            bat_parts_minus = GetCurrentParts(sol_lb,"I_minus",numBatParts) 
#        print ("Bat parts plus:",bat_parts_plus)
#        print ("Bat parts minus:",bat_parts_minus)
        #lb += added_obj
#        print ('iter',k, "lb:",lb)
        if lb > incumbent:
            incumbent = lb*1.0
            gap = (obj-incumbent)/incumbent
#            print ('updated gap:', gap)
            elapsed = time.time() - clock
            iterfile = open("iterations-method-"+method+"-scen"+scenario+"-m"+str(numSocParts)+"-n"+str(numBatParts)+".csv",'a')
            iterfile.write(str(elapsed)+","+str(incumbent)+","+str(obj)+","+str(gap)+"\n")
#            iterfile.close()  
        if k % ubCheck == 0:
            mean_soc = scipy.average([s["B_soc_end"] for s in sol_lb])
            mean_soc = 0.5
#            print ('mean_soc',mean_soc)
            solschecked=0
            for sol in sol_lb: 
                elapsed = time.time()-clock
                if elapsed > timeLimit: break
                if solschecked == solsToCheck: break
                for idx in range(len(sol["X"])):
                    sol["X"][idx] = (
                        sol["X"][idx][0],
                        float(int(sol["X"][idx][1]+2.5)/5)*5
                        ) #round to nearest 5 panels
                if str(sol["W"])+str(sol["X"]) not in ub_sols:
                    solschecked += 1
#                    print ("attempting new candidate.", str(sol["W"])+str(sol["X"]))
                    Partition_Univariate.OutputAllSols([sol],output_dir = IFS_dir, scenario_name = scenario)
                    if method == "nagarajan":
                        objU,solU,frac_over = Partition_Nagarajan.RunParallelUBDay(scenario,num_days=num_days,probSize=probSize,batPartsPlus=bat_parts_plus,
                            batPartsMinus=bat_parts_minus, boundary_soc=mean_soc, socParts=socParts)     
                    elif method == "univariate":
                        objU,solU,frac_over = Partition_Univariate.RunParallelUBDay(scenario,num_days=num_days,probSize=probSize,batPartsPlus=bat_parts_plus,
                            batPartsMinus=bat_parts_minus, boundary_soc=mean_soc) 
                    elif method == "gounaris":
                        objU,solU,frac_over = Partition_Univariate_G.RunParallelUBDay(scenario,num_days=num_days,probSize=probSize,batPartsPlus=bat_parts_plus,
                            batPartsMinus=bat_parts_minus, boundary_soc=mean_soc) 
                    elif method == "bivariate":
                        objU,solU,frac_over = Partition_Bivariate.RunParallelUBDay(scenario,num_days=num_days,probSize=probSize,batPartsPlus=bat_parts_plus,
                            batPartsMinus=bat_parts_minus, boundary_soc=mean_soc, socParts=socParts) 
#                    print ("SOL MIP UB",objU)
                    ub_sols.append(str(sol["W"])+str(sol["X"]))
                    frac_overs.append(frac_over)
                    if objU < obj:
                        if MINLP:
                            Partition_Univariate.OutputAllSols(solU,IFS_dir,scenario)
                            MINLP = MIPtoMINLP.MINLPSolution(inputs_dir,IFS_dir,scenario,time_horizon=num_days*probSize)
                            try: 
                                MINLP.RunSolution(num_intervals=num_days)
                                objU = MINLP.CalculateTotalCost()
                            except AttributeError:
                                objU *= 1.0
                            except AssertionError:
#                                print ("error in code. aborted solution.")
                                objU = scipy.infty
                    if objU <  obj:
#                        print ("new incumbent upper bound found. obj:",objU)
                        obj = objU
                        sol_ub = solU
                        gap = (obj-incumbent)/incumbent
#                        print ('updated gap:', gap)
                        elapsed = time.time() - clock
                        iterfile.write(str(elapsed)+","+str(incumbent)+","+str(obj)+","+str(gap)+"\n")
#                        iterfile.close()  
                        if gap < tol:
                            return incumbent,obj,gap,(k),scipy.average(frac_overs)
#            if sub != str(sol_ub) and gap >= tol:
#                sol_ub, obj, gap = FindBestSOC(sol_ub, obj, gap, incumbent, 
#                    scenario,num_days=num_days,
#                    probSize=probSize,batPartsPlus=bat_parts_plus,
#                    batPartsMinus=bat_parts_minus, MINLP=MINLP)    
        iterfile.close()
        k += 1
        lambs = UpdateLambdaPH(lambs,sol_lb,rho_bat)
        design_mults, design_objs, added_obj = UpdateDesignMults(design_mults,
                sol_lb, rho, proximal_term=False)
    return incumbent,obj,gap,(k-1),scipy.average(frac_overs)
    
def FindBestSOC(sol_ub, obj, gap, incumbent, scenario,num_days=365,
        probSize=24,batPartsPlus=[0,1], batPartsMinus=[0,1],
        tol=0.05, MINLP=False):
    """
    uses bisection on battery SOC at boundary points to obtain new optimal 
    solution.
    """
    soc_mid = 0.5
    Partition_Univariate.OutputAllSols([sol_ub[0]],output_dir = "./../sol_files/"+scenario+"/", scenario_name = scenario)
    o_mid,sol_mid,frac_over = Partition_Univariate.RunParallelUBDay(scenario,
        num_days=num_days,probSize=probSize,batPartsPlus=batPartsPlus,
        batPartsMinus=batPartsMinus, boundary_soc=soc_mid)
    dist = 0.25
    while dist > tol:
        soc_hi = soc_mid + dist
        soc_lo = soc_mid - dist
        o_hi,s_hi,frac_over = Partition_Univariate.RunParallelUBDay(scenario,
            num_days=num_days,probSize=probSize,batPartsPlus=batPartsPlus,
            batPartsMinus=batPartsMinus, boundary_soc=soc_hi)
        o_lo,s_lo,frac_over = Partition_Univariate.RunParallelUBDay(scenario,
            num_days=num_days,probSize=probSize,batPartsPlus=batPartsPlus,
            batPartsMinus=batPartsMinus, boundary_soc=soc_lo)
        if o_lo < o_hi and o_lo < o_mid:
            soc_mid = soc_lo
            o_mid = o_lo
            sol_mid = s_lo
        elif o_hi < o_mid:
            soc_mid = soc_hi
            o_mid = o_hi
            sol_mid = s_hi
        dist /= 2
    if MINLP:
        Partition_Univariate.OutputAllSols(sol_mid,output_dir = scenario+"/", scenario_name = scenario)
        MINLP_sol = MIPtoMINLP.MINLPSolution(
            "./../inputs/",
            "./../sol_files/"+scenario+"/",
            scenario,
            time_horizon=num_days*probSize)
        try: 
            MINLP_sol.RunSolution(num_intervals=num_days,output=True)
            o_mid = MINLP.CalculateTotalCost()
        except AttributeError:
            print ("Infeasible UB solution.")
            o_mid = scipy.infty
        except AssertionError:
            print ("error in code. aborted solution.")
            o_mid = scipy.infty
    if o_mid < obj:
        sol_ub = sol_mid
        obj = o_mid
        gap = (obj-incumbent)/incumbent
    return sol_ub, obj, gap
        
        
    
    
def GetMostPopularDesign(sols):
    """ Returns the most popular design decision when given a collection of
    solutions to the lower bound subproblems as input.
    sols -- solutions to lower bound subproblems
    retval -- a tuple of the binary/integer solutions and the top frequency"""
    solutions = {}
    Ws = []
    Xs = []
    for sol in sols:
        if str(sol["W"])+str(sol["X"]) in solutions.keys():
            solutions[str(sol["W"])+str(sol["X"])] += sol["end"]-sol["start"]+1
        else:
            solutions[str(sol["W"])+str(sol["X"])] = sol["end"]-sol["start"]+1
            Ws.append(sol["W"])
            Xs.append(sol["X"])
    incumbent = None
    hiFreq = 0
    for i in range(len(Ws)):
        if solutions[str(Ws[i])+str(Xs[i])] > hiFreq:
            hiFreq = solutions[str(Ws[i])+str(Xs[i])]
            incumbent = (Ws[i],Xs[i])
    return incumbent, hiFreq
        
def UpdateLambda(lambs,sol_lb,ubs,t):
    sol = sorted((x["start"],x) for x in sol_lb)
    soc_starts = scipy.array([x[1]["B_soc_start"] for x in sol[1:]])
    soc_ends = scipy.array([x[1]["B_soc_end"] for x in sol[:-1]])
    lbs = scipy.array([x[1]["lb"] for x in sol])
    #print 'lbs',lbs
    #print 'socStart', soc_starts
    #print 'socEnd', soc_ends
    #calculate step sizes
    stepsizes = t * (ubs[1:]+ubs[:-1]-lbs[1:]-lbs[:-1])
    new_lambs = scipy.maximum(scipy.zeros_like(soc_starts),scipy.array(lambs[1:-1]) + stepsizes*(soc_starts-soc_ends))
    return [0.0] + list(new_lambs) + [0]

def UpdateLambdaPH(lambs,sols,multiplier):
    sols = sorted((x["start"],x) for x in sols)
    boundary_invs = scipy.array([x[1]["battery_inv"] for x in sols])
    boundary_inv_mean = scipy.average(boundary_invs)
    lambs += multiplier * (boundary_invs - boundary_inv_mean)
    return lambs
        
def UpdateDesignMults(design_mults,sols,rho,proximal_term=False):
    """
    Updates design decision Lagrange multipliers based on the solution given.
    design_mults -- current design decision multipliers (for iteration k)
    sols -- solutions from iteration k of lower bound problem
    rho -- scalars used to multiply by difference between a single sol and 
            the mean
    num_periods -- number of time periods in full problem, used to  
            determine the mean purchase decision, amortized over time
    proximal -- include proximal term in objective function
    
    retval -- new design decision multipliers (for iteration k+1), objective
            values including proximal term, and added objective value due to 
            proximal term
    """
    #determine the mean occurence of each part of the design.
    num_periods = sum([sol["end"] - sol["start"] + 1 for sol in sols])
    means = [scipy.zeros_like(design_mults[0][0]),scipy.zeros_like(design_mults[0][1])]
    for s in range(len(sols)):
        for i in range(len(sols[s]["W"])):
            means[0][i] += 1.0*sols[s]["W"][i][1]*(sols[s]["end"]-sols[s]["start"]+1)/num_periods
        for i in range(len(sols[s]["X"])):
            #print sols[s]["X"][i]
            means[1][i] += 1.0*sols[s]["X"][i][1]*(sols[s]["end"]-sols[s]["start"]+1)/num_periods
    #print sols[0]["W"], sols[0]["X"]
    #print "means:"
    #for mean in means: print mean,"\n\n"
    #now determine the weight shift.
    #print "old design mults:", design_mults
    for s in range(len(sols)):
        for i in range(len(design_mults[s][0])):  #W variables (generators, batteries)
            design_mults[s][0][i] += rho[0][i] * (sols[s]["W"][i][1]-means[0][i])
        for i in range(len(design_mults[s][1])):  #X variables (pv systems)  
            design_mults[s][1][i] += rho[1][i] * (sols[s]["X"][i][1]-means[1][i])
    design_objs = [(x*1.0,y*1.0) for x,y in design_mults]
    added_obj = 0.0
    if proximal_term:
        for s in range(len(sols)):
            for i in range(len(design_mults[s][0])):
                added_obj += (rho[0][i]/2.0) * means[0][i]**2
                design_objs[s][1][i] += (rho[0][i]/2.0) * ((1-means[0][i])**2 - means[0][i]**2) 
            for i in range(len(design_mults[s][1])):
                added_obj += (rho[1][i]/2.0) * means[1][i]**2
                design_objs[s][1][i] += (rho[1][i]/2.0) * ((1-means[1][i])**2 - means[1][i]**2)   
        #print s,"new design mults:", design_mults 
    #for s in design_mults: 
    #    print "mults:", s
    return design_mults, design_objs, added_obj
    
def GenerateDesignMultipliers(scenario="ll1",
                tech_filename = "./../inputs/Opt_Model_J_Inputs.csv",
                max_filename = "./../inputs/Opt_Model_J_MAX.csv",
                pv_spec_filename = "./../inputs/Opt_Model_PV_Inputs.csv",
                num_subproblems = 365,
                probSize = 24,
                proximal = False):
    """Creates an initial set of Lagrange multipliers for the decision variables
    in the lower bound problem.
    scenario -- location indicator
    tech_filename -- file location for technology inputs
    max_filename -- file location for limits of number of each technology
    pv_spec_filename -- file location for pv system specs 
    num_subproblems -- number of subproblems to solve
    
    retval -- a zero for every design decision variable, for every problem, in
    a list of 2-tuples of arrays (one for W, one for X).
    """
    #imports and method definitions
    import InputReader
    def W_var(tech,k): return "W_"+tech.GetName()+".K"+str(k)
    def X_var(pv): return "X_"+pv.GetName() 
    def flatten(l):
        out = []
        for item in l:
            if isinstance(item, (list, tuple)):
                out.extend(flatten(item))
            else:
                out.append(item)
        return out
    
    #weight to apply to each purchase decision cost increment.
    weight = (1.0*probSize)/ (8760.0 * 3.0)
    weight_battery = (1.0*probSize)/ (8760.0 * 3.0) 
    
    #load inputs in a form so that we can make variables from them
    technologies = InputReader.ReadTechnologies(tech_filename, 
                                0.8,8760)
    tech_maxes = InputReader.ReadTechMaxes(max_filename)
    pv_specs = InputReader.ReadPVSystemSpecs(pv_spec_filename)
    #print pv_specs
    for tech in technologies:
        if tech.GetType() != "PVArray":
            tech.SetTechMax(tech_maxes[scenario][tech.GetName()])
    W_names = flatten([["W_" + tech.GetName() + ".K" + str(i+1) 
            for i in range(tech_maxes[scenario][tech.GetName()])] 
            for tech in technologies 
                if tech.GetType() != "PVArray"])
    W_mults = scipy.array(flatten([[tech.GetPurchaseCost()
            for i in range(tech_maxes[scenario][tech.GetName()])] 
            for tech in technologies 
                if tech.GetType() != "PVArray"])) * weight
    max_battery_cost = max(flatten([[tech.GetPurchaseCost()
            for i in range(tech_maxes[scenario][tech.GetName()])] 
            for tech in technologies 
            if tech.GetType() == "Battery"]))
    max_battery_cap = max(flatten([[tech.GetMaxPower()
            for i in range(tech_maxes[scenario][tech.GetName()])] 
            for tech in technologies 
            if tech.GetType() == "Battery"]))
    #print "max bat cap:", max_battery_cap
    X_names = ["X_" + pv_specs[pv]["name"] for pv in pv_specs.keys() 
            if pv_specs[pv]["include"] > 0.5]
    X_mults = scipy.array([pv_specs[pv]["purchase_cost"] for pv in pv_specs.keys()]) * weight
    W_vals = scipy.zeros(len(W_names))
    X_vals = scipy.zeros(len(X_names))
    #print "RHO_W:",W_mults
    #print "RHO_X:",X_mults 
    return [(W_vals,X_vals) for i in range(num_subproblems)], (W_mults,X_mults), max_battery_cost*weight_battery/max_battery_cap
    

def GetSolVal(name,l):
    for pair in l:
        if pair[0] == name: return pair[1]
    return None
      
def GetLinearSol(scenario):
    import InputReader
    direc = "./../inputs/"
    IFS_dir = "./../sol_files/"+scenario+"/"
    IFS = InputReader.GetIFS(IFS_dir)
    scalar_filename = direc + "Opt_Model_Scalars.csv"
    load_filename = direc + "Opt_Model_Demand.csv"
    fc_filename = direc + "Opt_Model_Fuel.csv"
    tech_filename = direc + "Opt_Model_J_inputs.csv"
    pv_filename = direc + "Opt_Model_Solar.csv"
    pv_spec_filename = direc + "Opt_Model_PV_inputs.csv"
    max_filename = direc + "Opt_Model_J_MAX.csv"  
    cuts_filename = direc + "Opt_Model_Cuts.csv" 
    mincap_filename = direc + "mingencap.csv" 
    IFS = InputReader.GetIFSVals(IFS_dir)
    inputs = InputReader.GetInputs(scalar_filename,load_filename,
            pv_filename,pv_spec_filename,fc_filename,tech_filename,
            max_filename,cuts_filename,8760,mincap_filename=mincap_filename) 
    InputReader.setBounds(inputs,scenario)
    sol = Partition_Univariate.SolvePartition(inputs, scenario, 1, 8760, 
        inv_mult=0.0, upper_bound = False, mipgap=0.05, design_mults = [0,0],
        batPartsPlus = [0.0,1.0], batPartsMinus=[0.0,1.0],
        reslim=3600, IFS = IFS, output = True, 
        min_PV = 0, amortize_on_demand = False, scen=1, nscens=1)
    gap = (sol["obj"]-sol["lb"])/sol["lb"]
    return sol["lb"], sol["obj"], gap, 0
            
if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1: 
        scenario = sys.argv[1]
        print ("Scenario selected: ", scenario)
    else: 
        print ("No scenario argument provided; using ll1.")
        scenario = "ll1"
    if len(sys.argv) > 3: 
        numBatParts = int(sys.argv[3])
        print ("Current partitions selected: ", numBatParts)
    else: 
        print ("No argument for current partitions provided; using default of 4.")
        numBatParts = 4
    if len(sys.argv) > 2: numSocParts = int(sys.argv[2])
    else: 
        print ("No argument for SOC partitions provided; using default of 1.")
        numSocParts = 1
    if len(sys.argv) > 4: 
        method = sys.argv[4]
        print ("partitioning model selected: ", method)
    else: 
        print ("No argument for SOC partitions provided; using 'univariate' or model (U)")
        method = "univariate"
    if len(sys.argv) > 5: 
        if sys.argv[5] == 'minlp':
            MINLP = True
        else: 
            print ("argument for MINLP solve not recognized; solving MILP only.")
            MINLP = False
    else: 
        print ("no argumemt provided for MINLP solve; solving MILP only.")
    print ("scenario:",scenario)
    tech_filename = "./../inputs/Opt_Model_J_inputs.csv"
    max_filename = "./../inputs/Opt_Model_J_MAX.csv"
    pv_spec_filename = "./../inputs/Opt_Model_PV_Inputs.csv"
    outfile = open("results_method-"+method+"_m"+str(numSocParts)+"_n"+str(numBatParts)+"_"+scenario+".csv","w") 
    start = time.time()
    lb,obj,gap,iters,frac_over = RunPH(scenario,num_days=365,probSize=24,
            itermax=1000,tol=0.05,n=500,ubCheck=1,
            tech_filename = tech_filename,
            max_filename = max_filename,
            pv_spec_filename = pv_spec_filename,
            lbFind=1, numSocParts=numSocParts, numBatParts=numBatParts, MINLP=False,nonuniform_part=False,
            findgen=True, timeLimit=900, solsToCheck=5, method=method
            )
    #
    elapsed = (time.time() - start) #measure time elapsed in minutes
    outfile.write("Scenario,LB,UB,gap,iters,time,fracOverLim\n")
    outfile.write(scenario+","+str(lb)+","+str(obj)+","+str(gap)+","+str(iters)+","+str(elapsed)+","+str(frac_over)+","+"\n")
    print("\n\nTerminated. Results:")
    print("method: ",method)
    print("scenario: ",scenario)
    print("lb: ",lb)
    print("obj: ",obj)
    print("gap: ",gap)
    print("total time:",elapsed)
    print("\n\n\n")
    
   