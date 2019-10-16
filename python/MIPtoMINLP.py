"""
This module takes in as input a set of technologies (via InputReader) and
a MILP-feasible solution and outputs a MINLP-feasible solution.  This is done 
as an attempt to get a valid upper bound on the nonlinear model.
"""

import InputReader
import scipy

class MINLPSolution(object):
    def __init__(self,inputs_dir,IFS_dir,scenario,time_horizon=8760):
        self.scenario = scenario
        self.IFS_dir = IFS_dir
        self.inputs = self.GetAllInputs(inputs_dir, time_horizon)
        self.time_horizon = time_horizon
        self.technologies = self.inputs["technologies"]
        #print self.technologies[-1]
        self.MIPSol = self.ImportMIPSolution(IFS_dir)
        self.StartSOC = self.inputs["scalars"]["b_init"]
        self.StartSOC = 0.5
        self.spinning_reserve_factor = self.inputs["scalars"]["k_sol"]
        self.fuel_costs = scipy.array(self.inputs["fuel_costs"][self.scenario])
        self.generators = self.GetGenerators()
        try:
            self.battery = self.GetBattery()
        except KeyError:
            self.battery = None
        self.pv_array = self.GetPVArray()
        self.fuel_as = scipy.array([gen.GetFuelUseA() for gen in self.generators]) / 1000000
        self.fuel_bs = scipy.array([gen.GetFuelUseB() for gen in self.generators]) / 1000
        self.fuel_cs = scipy.array([gen.GetFuelUseC() for gen in self.generators])

        self.gen_caps = scipy.array([gen.GetMaxPower() for gen in self.generators])
        self.tau = self.inputs["scalars"]["tau"]
        self.nu = self.inputs["scalars"]["nu"]
        self.totalCost = self.MIPSol["Z"][self.scenario]
        #Initialize currents
        self.B_plus = 0
        self.B_minus = 0
        self.I_plus = 0.0
        self.Z_plus = 0.0
        self.I_minus = 0.0
        self.Z_minus = 0.0
        self.P_solar = 0.0
        self.reserve = 0.0
        self.threshold = 0.0
        self.F_tilde = 0.0
        #Storing solution in arrays.
        self.Gs = scipy.array([[0 for g in self.generators] for t in range(self.time_horizon)])
        self.P_out_gs = scipy.array([[0.0 for g in self.generators] for t in range(self.time_horizon)])
        self.P_out_bs = scipy.array([0.0 for t in range(self.time_horizon)])
        self.P_in_bs = scipy.array([0.0 for t in range(self.time_horizon)])
        self.P_solars = scipy.array([0.0 for t in range(self.time_horizon)])
        self.B_pluses = scipy.array([0.0 for t in range(self.time_horizon)])
        self.B_minuses = scipy.array([0.0 for t in range(self.time_horizon)])
        self.I_pluses = scipy.array([0.0 for t in range(self.time_horizon)])
        self.I_minuses = scipy.array([0.0 for t in range(self.time_horizon)])
        self.Z_pluses = scipy.array([0.0 for t in range(self.time_horizon)])
        self.Z_minuses = scipy.array([0.0 for t in range(self.time_horizon)])
        self.B_socs = scipy.array([0.0 for t in range(self.time_horizon)])
        self.F_tildes = scipy.array([0.0 for t in range(self.time_horizon)])
        self.eps = 1e-8
        
    def Reset(self):
        self.B_plus = 0
        self.B_minus = 0
        self.I_plus = 0.0
        self.Z_plus = 0.0
        self.I_minus = 0.0
        self.Z_minus = 0.0
        self.P_solar = 0.0
        self.reserve = 0.0
        self.threshold = 0.0
        self.F_tilde = 0.0
        self.Gs = scipy.array([[0 for g in self.generators] for t in range(self.time_horizon)])
        self.P_out_gs = scipy.array([[0.0 for g in self.generators] for t in range(self.time_horizon)])
        self.P_out_bs = scipy.array([0.0 for t in range(self.time_horizon)])
        self.P_in_bs = scipy.array([0.0 for t in range(self.time_horizon)])
        self.P_solars = scipy.array([0.0 for t in range(self.time_horizon)])
        self.B_pluses = scipy.array([0.0 for t in range(self.time_horizon)])
        self.B_minuses = scipy.array([0.0 for t in range(self.time_horizon)])
        self.I_pluses = scipy.array([0.0 for t in range(self.time_horizon)])
        self.I_minuses = scipy.array([0.0 for t in range(self.time_horizon)])
        self.Z_pluses = scipy.array([0.0 for t in range(self.time_horizon)])
        self.Z_minuses = scipy.array([0.0 for t in range(self.time_horizon)])
        self.B_socs = scipy.array([0.0 for t in range(self.time_horizon)])
        self.F_tildes = scipy.array([0.0 for t in range(self.time_horizon)])
        
    def CalculateTotalCost(self,output=False):
        if output:
            print "Purchase", self.CalculatePurchaseCost()
            print "LifeCycles", self.CalculateLifecycleCost()
            print "Fuel", self.CalculateFuelCost()
            print "Weight", self.CalculateWeightCost() 
            print "volume", self.CalculateVolumeCost()
            print "Total", (self.CalculatePurchaseCost() + 
                    self.CalculateLifecycleCost() + 
                    self.CalculateFuelCost() + self.CalculateWeightCost() + 
                    self.CalculateVolumeCost())
        return (self.CalculatePurchaseCost() + self.CalculateLifecycleCost() + 
                self.CalculateFuelCost() + self.CalculateWeightCost() + 
                self.CalculateVolumeCost())
        
    def CalculatePurchaseCost(self):
        return (sum([gen.GetPurchaseCost() for gen in self.generators]) + 
                self.battery.GetPurchaseCost() + 
                self.pv_array.GetPurchaseCost())*self.time_horizon/8760.0     
    
    def CalculateFuelUse(self):
        self.F_tilde = scipy.sum(self.fuel_as * self.P_out_g * self.P_out_g 
                + self.fuel_bs * self.P_out_g + self.fuel_cs * self.G, 1)
                
    def GetBatteryLifecycles(self):
        return ( ( scipy.sum(self.I_pluses)*self.battery.GetASOC() - 
                scipy.sum(self.Z_pluses)*self.battery.GetDSOC() ) / 
                self.battery.GetReferenceCapacity() )
                                                                                                
    def CalculateFuelCost(self):
        return scipy.sum(self.fuel_costs * self.F_tildes)
        
    def CalculateLifecycleCost(self):
        gen_cycle_costs = scipy.array([gen.GetLifeCycleCost() for gen in self.generators])
        bat_cycles = self.GetBatteryLifecycles()
        return ( scipy.sum(gen_cycle_costs * self.Gs) + 
                self.battery.GetLifeCycleCost()*bat_cycles)
        
    def CalculateWeightCost(self): return 0.0
    
    def CalculateVolumeCost(self): return 0.0
    
    def SetThreshold(self,threshold):
        self.threshold = threshold
                        
    def GetAllInputs(self,inputs_dir,time_horizon):
        """Gets all inputs to the problem, including technology specs, fuel 
        costs, etc."""
        scalar_filename = inputs_dir+"Opt_Model_Scalars.csv"
        load_filename = inputs_dir+"Opt_Model_Demand.csv"
        fc_filename = inputs_dir+"Opt_Model_Fuel.csv"
        tech_filename = inputs_dir+"Opt_Model_J_inputs.csv"
        pv_filename = inputs_dir+"Opt_Model_Solar.csv"
        pv_spec_filename = inputs_dir+"Opt_Model_PV_Inputs.csv"
        max_filename = inputs_dir+"Opt_Model_J_MAX.csv"  
        cuts_filename =  inputs_dir+"Opt_Model_Cuts.csv"
        #read and initialize inputs
        inputs = InputReader.GetInputs(scalar_filename,load_filename,pv_filename,
                pv_spec_filename,fc_filename,tech_filename,max_filename,
                cuts_filename,time_horizon)
        return inputs
        
        
    def ImportMIPSolution(self,IFS_dir):
        """Retrieves as a dictionary the solution retrieved by solving the MIP
        relaxation of the model in CPLEX, by reading the set of files generated
        by the output.
        IFS_dir - directory containing IFS
        retval - dictionary with all nonzero values
        """
        solution = {}
        solution["W"] = InputReader.ReadInitialSols3Index(IFS_dir+"WW.csv")
        solution["X"] = InputReader.ReadInitialSols2Index(IFS_dir+"X.csv")
        solution["Z"] = InputReader.ReadInitialSols1Index(IFS_dir+"Z.csv")
        solution["B_soc"] = InputReader.ReadInitialSols4Index(IFS_dir+"B_soc.csv")
        solution["P_out"] = InputReader.ReadInitialSols4Index(IFS_dir+"P_out.csv")
        solution["P_solar"] = InputReader.ReadInitialSols3Index(IFS_dir+"P_solar.csv",style='f')
        solution["P_in"] = InputReader.ReadInitialSols4Index(IFS_dir+"P_in.csv")
        solution["B_plus"] = InputReader.ReadInitialSols4Index(IFS_dir+"B_plus.csv",style='i')
        solution["B_minus"] = InputReader.ReadInitialSols4Index(IFS_dir+"B_minus.csv",style='i')
        solution["I_plus"] = InputReader.ReadInitialSols4Index(IFS_dir+"I_plus.csv")
        solution["I_minus"] = InputReader.ReadInitialSols4Index(IFS_dir+"I_minus.csv")
        solution["Z_plus"] = InputReader.ReadInitialSols4Index(IFS_dir+"Z_plus.csv")
        solution["Z_minus"] = InputReader.ReadInitialSols4Index(IFS_dir+"Z_minus.csv")
        solution["F_tilde"] = InputReader.ReadInitialSols2Index(IFS_dir+"F_tilde.csv") 
        solution["G"] = InputReader.ReadInitialSols4Index(IFS_dir+"G.csv",style="i") 
        #print solution["G"]
        return solution
    
    def GetGenerators(self):
        gens = []
        for techname in self.MIPSol["G"][self.scenario].keys():
            if techname[0] == "G":
                i=1
                for twin in self.MIPSol["G"][self.scenario][techname].keys():
                    gens.append(self.GetTech(techname).CopyDieselGenerator(i))
                    i += 1
        sorted_gens = sorted([(gen.GetName(),gen.GetTwin(),gen) for gen in gens])
        return [x[2] for x in sorted_gens]
    
    def GetTech(self,name):
        for tech in self.technologies:
            if tech.GetName() == name:  return tech
        return -1
        
        
    def GetBattery(self):
        for tech in self.technologies:
            if tech.GetName() in self.MIPSol["B_plus"][self.scenario].keys() and tech.GetType() == "Battery":
                return tech.CopyBattery(1)
        return -1
    
    def GetPVArray(self):
        for tech in self.technologies:
            if tech.GetName() in self.MIPSol["X"][self.scenario].keys() and tech.GetType() == "PVArray":
                num_panels = self.MIPSol["X"][self.scenario][tech.GetName()]
                return tech.CopyPVArray(1,num_panels)    
        return -1
        
    def GetMIPSolutionValue(self,output,techname,twin,t):
        """Gets the MIP solution value at time period t for some defined output.
        output -- output desired (e.g., P_out)
        techname -- technology name (usin set notation, e.g., "G1","G2","B1" 
        twin -- twin number for generators or batteries, e.g., "K1", K2"
        t -- time period number, e.g., "TT1"
        retval -- floating point number indicating the value in the solution"""
        try:
            if output in ["G","B_plus","B_minus"]:
                #print output, self.MIPSol[output][self.scenario].keys()
                return int(self.MIPSol[output][self.scenario][techname][twin]["TT"+str(t)])
            else:
                return float(self.MIPSol[output][self.scenario][techname][twin]["TT"+str(t)])
        except KeyError, e:
            #print repr(e)
            return 0.0
            
    def GetMIPFuel(self,t):
        """Gets the fuel used according to the MIP solution at time period t."""
        try:
            return float(self.MIPSol["F_tilde"][self.scenario]["TT"+str(t)])
        except KeyError:
            return 0.0
            
    def GetMIPSolar(self,t):
        """Gets the PV power out according to the MIP solution at period t."""
        try:
            return float(self.MIPSol["P_solar"][self.scenario][self.pv_array.GetName()]["TT"+str(t)])
        except KeyError:
            return 0.0
            
    def GetIntervalInputs(self,start_period,end_period):
        """obtain inputs for a reset interval (zeros for current values)"""
        num_periods = end_period-start_period+1
        self.G = scipy.array([[self.GetMIPSolutionValue("G",gen.GetName(),"K"+str(gen.GetTwin()),t) for gen in self.generators] for t in range(start_period,end_period+1)])
        self.P_out_g = scipy.array([[self.GetMIPSolutionValue("P_out",gen.GetName(),"K"+str(gen.GetTwin()),t) for gen in self.generators] for t in range(start_period,end_period+1)])
        self.P_out_b = scipy.array([self.GetMIPSolutionValue("P_out",self.battery.GetName(),"K1",t) for t in range(start_period,end_period+1)])
        self.P_in_b = scipy.array([self.GetMIPSolutionValue("P_in",self.battery.GetName(),"K1",t) for t in range(start_period,end_period+1)])
        self.B_plus = scipy.array([self.GetMIPSolutionValue("B_plus",self.battery.GetName(),"K1",t) for t in range(start_period,end_period+1)])
        self.B_minus = scipy.array([self.GetMIPSolutionValue("B_minus",self.battery.GetName(),"K1",t) for t in range(start_period,end_period+1)])
        self.P_solar = scipy.array([self.GetMIPSolar(t) for t in range(start_period,end_period+1)])
        self.I_minus = scipy.zeros(num_periods)
        self.I_plus = scipy.zeros(num_periods)
        self.Z_minus = scipy.zeros(num_periods)
        self.Z_plus = scipy.zeros(num_periods)
        self.reserve = self.P_solar * self.spinning_reserve_factor
        self.SOC_target = scipy.array([self.GetMIPSolutionValue("B_soc",self.battery.GetName(),"K1",t) for t in range(start_period,end_period+1)])
        self.SOC = self.SOC_target*1.0
        if start_period == 1:
            self.SOC0 = self.StartSOC
        else:
            self.SOC0 = self.SOC_target[-1] * 1.0
            
    def ComputeIntervalSOCs(self):
        for i in range(len(self.SOC)):
            self.SOC[i] = self.GetRealSOC(i)
        
    def RunDispatchPeriods(self,start_period,end_period):
        """
        Takes as input dispatch for full reset cycle, runs nonlinear model
        to obtain battery state of charge, assesses the required diesel 
        generation (and deadline period), and runs additional diesel generators
        as needed.
        """
        #print start_period
        self.GetIntervalInputs(start_period,end_period)
        self.ComputeIntervalSOCs()
        neg_bat = (self.SOC.min() < 0.0)
        hi_bat = (self.SOC.max() > 1.0)
        #if self.P_solar.max() > 0.001:
        #    first_solar = scipy.argmax(self.P_solar > 0)
        #else:
        #    first_solar = end_period-start_period+1
        #neg = 0
        while neg_bat or hi_bat or abs(self.SOC_target[-1] - self.SOC[-1]) > self.eps or not self.SpinningReserveIsAllMet():
            #print "maindisp"
            self.FixNegHiBats()
            neg_bat = (self.SOC.min() < -1*self.eps)
            hi_bat = (self.SOC.max() > 1.0+self.eps)
            #print "NH:", self.SOC#self.B_minus[-1], self.B_plus[-1], self.P_out_b[-1], self.P_in_b[-1], self.P_out_g[-1]
            self.CheckSpinningReserve()
            #print "SPIN:", self.SOC#self.B_minus[-1], self.B_plus[-1], self.P_out_b[-1], self.P_in_b[-1], self.P_out_g[-1]
            #print "FLOWS:", self.SOC#self.B_minus[-1], self.B_plus[-1], self.P_out_b[-1], self.P_in_b[-1], self.P_out_g[-1]
            self.CheckBatteryFlows()
            self.ComputeIntervalSOCs()
            #print "COMPUTE:", self.B_minus[-1], self.B_plus[-1], self.P_out_b[-1], self.P_in_b[-1], self.P_out_g[-1]
            self.MatchTargetSOC(-1,self.SOC_target[-1])
            self.CheckBatteryFlows()
            self.ComputeIntervalSOCs()
            neg_bat = (self.SOC.min() < -1*self.eps)
            hi_bat = (self.SOC.max() > 1.0+self.eps)
            #print start_period, neg_bat, hi_bat, self.SOC_target[-1], self.SOC[-1], self.SpinningReserveIsAllMet()
            #print "MATCH", self.SOC#self.B_minus[-1], self.B_plus[-1], self.P_out_b[-1], self.P_in_b[-1], self.P_out_g[-1]
        self.CalculateFuelUse()
        self.StoreDispatch(start_period,end_period)
    
    def CheckBatteryFlows(self):
        redo = True
        while redo:
            redo = False
            #print "flow check"
            for idx in range(len(self.G)):
                if idx == 0:
                    self.battery.SetSOC(self.SOC0)
                else: 
                    self.battery.SetSOC(self.SOC[idx-1])
                if self.B_minus[idx] == 1:
                    while self.P_out_b[idx] > self.battery.GetMaxDelivery(self.tau,self.nu)/self.battery.GetEfficiencyOut()+self.eps:
                        #print "flowoutoff", idx,self.SOC[idx],self.battery.GetMaxDelivery(self.tau,self.nu)/self.battery.GetEfficiencyOut(),self.battery.GetSOC(),self.P_out_b[idx]
                        if idx == 0 or sum(self.G[idx]*self.gen_caps-self.P_out_g[idx]) >= self.P_out_b[idx] - self.battery.GetMaxDelivery(self.tau,self.nu)/self.battery.GetEfficiencyOut():
                            diesel_to_add = self.P_out_b[idx] - self.battery.GetMaxDelivery(self.tau,self.nu)
                            self.P_out_b[idx] = self.battery.GetMaxDelivery(self.tau,self.nu)/self.battery.GetEfficiencyOut()-self.eps
                            self.AddDieselOutput(diesel_to_add,idx)
                        else: 
                            over_discharge = self.P_out_b[idx] - self.battery.GetMaxDelivery(self.tau,self.nu)/self.battery.GetEfficiencyOut()
                            if over_discharge < self.RemainingDieselAvailable(idx):
                                self.P_out_b[idx] = self.battery.GetMaxDelivery(self.tau,self.nu)/self.battery.GetEfficiencyOut()
                                self.AddDieselOutput(over_discharge * self.battery.GetEfficiencyOut(), idx)
                            else:
                                P_deficit = over_discharge - self.RemainingDieselAvailable(idx)
                                I_deficit = P_deficit / self.battery.CalculateDischargeVoltage(self.tau)
                                new_soc = self.SOC[idx-1] + (I_deficit)/self.battery.GetReferenceCapacity()
                                assert new_soc < 1.0+self.eps and new_soc > -self.eps, "BAD SOC FLOWS "+str(new_soc)
                                self.MatchTargetSOC(idx-1,new_soc+self.eps)
                                redo = True
                        self.ComputeIntervalSOCs()
                if self.B_plus[idx] == 1:
                    while self.P_in_b[idx] > self.battery.GetMaxCharge(self.tau,self.nu)+self.eps:
                        #print "flowinoff", idx,self.SOC[idx],self.battery.GetMaxCharge(self.tau,self.nu),self.battery.GetSOC(),self.P_in_b[idx] 
                        if sum(self.P_out_g[idx]) >= self.P_in_b[idx] - self.battery.GetMaxCharge(self.tau,self.nu):
                            diesel_reduction = self.P_in_b[idx] - self.battery.GetMaxCharge(self.tau,self.nu)
                            self.P_in_b[idx] = self.battery.GetMaxCharge(self.tau,self.nu)-self.eps
                            self.ReduceDieselOutput(diesel_reduction,idx)
                        else:
                            self.ReduceDieselOutput(sum(self.P_out_g[idx]),idx)
                            self.P_in_b[idx] = self.battery.GetMaxCharge(self.tau,self.nu)
                        self.ComputeIntervalSOCs()
        
    def FixNegHiBats(self):
        neg_bat = (self.SOC.min() < -1*self.eps)
        hi_bat = (self.SOC.max() > 1.0+self.eps)
        while neg_bat or hi_bat:
            #print 'neghi'
            #print "neghi",neg_bat,hi_bat
            if neg_bat:
                first_neg_bat = scipy.argmax(self.SOC < 0)
                #print "first_neg_bat",first_neg_bat,self.SOC[scipy.argmax(self.SOC < 0)]
                self.MatchTargetSOC(first_neg_bat,self.eps)
                #print self.SOCF
            if hi_bat:
                first_hi_bat = scipy.argmax(self.SOC > 1.0)
                #print "first_hi_bat",first_hi_bat,self.SOC[scipy.argmax(self.SOC > 1.0)]
                self.MatchTargetSOC(first_hi_bat,1.0)
                #print self.SOC
            neg_bat = (self.SOC.min() < -1*self.eps)
            hi_bat = (self.SOC.max() > 1.0+self.eps)
        #print "NH", neg_bat, hi_bat
    
    def SpinningReserveIsAllMet(self):
        return all([self.SpinningReserveIsMet(self.SOC[idx],idx) for idx in range(len(self.G))])        
                            
    def CheckSpinningReserve(self):
        #add spinning reserve as needed.
        spinning_reserve_met = all([self.SpinningReserveIsMet(self.SOC[idx],idx) for idx in range(len(self.G))])
        tries = 1
        while not spinning_reserve_met:
            #print 'spin'
            for idx in range(len(self.G)):
                if not self.SpinningReserveIsMet(self.SOC[idx],idx):
                    deficit = self.reserve[idx] - self.GetSpinningReserve(self.SOC[idx],idx) 
                    if tries == 1:  
                         new_soc = self.SOC[idx-1] + (deficit)/self.battery.GetMaxPower()
                         self.MatchTargetSOC(idx-1,new_soc)
                    else:
                        self.AddDieselCapacityOnly(deficit,idx)
                    self.ComputeIntervalSOCs()
            tries += 1
            if tries % 10 == 0: print tries
            spinning_reserve_met = all([self.SpinningReserveIsMet(self.SOC[idx],idx) for idx in range(len(self.G))])           
    
    def MatchTargetSOC(self,idx,soc_target):
        #diesel_avail = scipy.sum(self.gen_caps - self.P_out_g[idx])
        #diesel_running = scipy.sum(self.P_out_g[idx])
        #print idx,soc_target
        if idx == 0:
            self.battery.SetSOC(self.SOC0) 
        else:
            self.battery.SetSOC(self.SOC[idx-1]) 
        if self.SOC[idx] > soc_target:
            diesel_reduction, delta_P_plus, delta_P_minus = self.GetDischargeDeficit(idx,soc_target)
            self.P_in_b[idx] += delta_P_plus
            self.P_out_b[idx] += delta_P_minus
            self.ReduceDieselOutput(diesel_reduction,idx)
        elif self.SOC[idx] < soc_target:
            diesel_to_add, delta_P_plus, delta_P_minus = self.GetChargeDeficit(idx,soc_target)
            self.P_in_b[idx] += delta_P_plus
            self.P_out_b[idx] += delta_P_minus
            self.AddDieselOutput(diesel_to_add,idx)
        self.ComputeIntervalSOCs()
        
            
    def GetDischargeDeficit(self,idx,soc_target):
        #print "dis", idx, soc_target, self.SOC[idx]
        if idx == 0:
            self.battery.SetSOC(self.SOC0) 
        else:
            self.battery.SetSOC(self.SOC[idx-1]) 
        if soc_target >= self.SOC[idx]: return 0,0,0
        if self.B_minus[idx] == 1:
            #if idx == len(self.G)-1 or idx == -1:
            #    print "DISDEF",self.SOC[idx],self.battery.GetMaxDelivery(self.tau,self.nu)/self.battery.GetEfficiencyOut(),self.battery.GetSOC(),self.P_out_b[idx]
            if self.SOC[idx-1] > soc_target:
                I_minus_increase = (self.SOC[idx]-soc_target)*self.battery.GetReferenceCapacity()#/self.battery.GetEfficiencyIn()
                P_minus_increase = I_minus_increase * self.battery.CalculateDischargeVoltage(self.tau)
                if P_minus_increase > self.battery.GetMaxDelivery(self.tau,self.nu)/self.battery.GetEfficiencyOut() - self.P_out_b[idx]: #"Discharge limit exceeded."
                    #print "max delivery exceeded."
                    assert idx != 0, "end of the line - no more discharge capacity."
                    P_deficit = P_minus_increase - (self.battery.GetMaxDelivery(self.tau,self.nu)/self.battery.GetEfficiencyOut() - self.P_out_b[idx])
                    #self.battery.SetSOC(self.SOC[idx-1])
                    if idx == 1:
                        self.battery.SetSOC(self.SOC0) 
                    else:
                        self.battery.SetSOC(self.SOC[idx-2])
                    I_deficit = P_deficit / self.battery.CalculateDischargeVoltage(self.tau)
                    new_target = self.SOC[idx-1] - I_deficit/self.battery.GetReferenceCapacity()
                    assert new_target <= 1.0+self.eps and new_target >= -1*self.eps, "invalid SOC target. prev soc: "+str(self.SOC[idx-1])+" target:"+str(new_target)
                    #print "next match d1:",idx-1,new_target-self.eps
                    self.MatchTargetSOC(idx-1,new_target-self.eps)
                    #print 'discharge target matched1:',new_target 
                    self.ComputeIntervalSOCs()
                    return self.GetDischargeDeficit(idx,soc_target)
                return P_minus_increase*self.battery.GetEfficiencyOut(), 0, P_minus_increase
            else:
                I_plus_increase = (soc_target-self.SOC[idx-1])*self.battery.GetReferenceCapacity()/self.battery.GetEfficiencyIn()
                P_plus_increase = I_plus_increase * self.battery.CalculateChargeVoltage(self.tau)
                if (P_plus_increase > self.battery.GetMaxCharge(self.tau,self.nu) - self.P_in_b[idx] or
                        P_plus_increase+self.P_out_b[idx]*self.battery.GetEfficiencyOut() > self.RemainingDieselAvailable(idx)):# "Charge limit exceeded."
                    assert idx != 0, "end of the line - no more charge capacity."
                    P_deficit = max(P_plus_increase - (self.battery.GetMaxCharge(self.tau,self.nu) - self.P_in_b[idx]),
                        P_plus_increase+self.P_out_b[idx]*self.battery.GetEfficiencyOut() - self.RemainingDieselAvailable(idx)) 
                    if idx == 1:
                        self.battery.SetSOC(self.SOC0) 
                    else:
                        self.battery.SetSOC(self.SOC[idx-2])
                    I_deficit = P_deficit / self.battery.CalculateChargeVoltage(self.tau)
                    new_target = self.SOC[idx-1] + self.battery.GetEfficiencyIn()*I_deficit/self.battery.GetReferenceCapacity()
                    assert new_target <= 1.0+self.eps and new_target >= -1*self.eps, "invalid SOC target. prev soc: "+str(self.SOC[idx-1])+" target:"+str(new_target)
                    #print "next match d2:",idx-1,new_target-self.eps
                    self.MatchTargetSOC(idx-1,new_target+self.eps)
                    #print 'discharge target matched2:',new_target, 'for period',idx-1 
                    self.ComputeIntervalSOCs()
                    return self.GetChargeDeficit(idx,soc_target)
                if P_plus_increase < 0:
                    #print "negative P_plus_increase."
                    return 0,0,0
                self.B_minus[idx] = 0
                self.B_plus[idx] = 1
                return P_plus_increase+self.P_out_b[idx]*self.battery.GetEfficiencyOut(), P_plus_increase, -1.0*self.P_out_b[idx] 
        elif self.B_plus[idx] == 1:   
            if self.SOC[idx-1] < soc_target:
                I_plus_reduction = (self.SOC[idx]-soc_target)*self.battery.GetReferenceCapacity()/self.battery.GetEfficiencyIn()
                P_plus_reduction  = I_plus_reduction * self.battery.CalculateChargeVoltage(self.tau)
                #assert P_plus_reduction+self.P_out_b[idx] < self.battery.GetMaxCharge(self.tau,self.nu): #"too much power in in period"+str(idx)
                return P_plus_reduction,-1.0*P_plus_reduction,0
            else:
                I_minus_increase = (self.SOC[idx-1]-soc_target)*self.battery.GetReferenceCapacity()
                P_minus_increase = I_minus_increase * self.battery.CalculateDischargeVoltage(self.tau)
                if P_minus_increase > self.battery.GetMaxDelivery(self.tau,self.nu)/self.battery.GetEfficiencyOut() - self.P_out_b[idx]: #"Discharge limit exceeded."
                    assert idx != 0, "end of the line - no more discharge capacity."
                    P_deficit = P_minus_increase - (self.battery.GetMaxDelivery(self.tau,self.nu)/self.battery.GetEfficiencyOut() - self.P_out_b[idx])
                    #self.battery.SetSOC(self.SOC[idx-1])
                    if idx == 1:
                        self.battery.SetSOC(self.SOC0) 
                    else:
                        self.battery.SetSOC(self.SOC[idx-2])
                    I_deficit = P_deficit / self.battery.CalculateDischargeVoltage(self.tau)
                    new_target = self.SOC[idx-1] - I_deficit/self.battery.GetReferenceCapacity()
                    assert new_target <= 1.0+self.eps and new_target >= -1*self.eps, "invalid SOC target. prev soc: "+str(self.SOC[idx-1])+" target:"+str(new_target)
                    #print "next match d3:",idx-1,new_target-self.eps
                    self.MatchTargetSOC(idx-1,new_target-self.eps)
                    #print 'discharge target matched3:',new_target 
                    self.ComputeIntervalSOCs()
                    return self.GetDischargeDeficit(idx,soc_target)
                self.B_plus[idx] = 0
                self.B_minus[idx] = 1
                return P_minus_increase*self.battery.GetEfficiencyOut()+self.P_in_b[idx], -1.0*self.P_in_b[idx], P_minus_increase  
        else:
            assert self.SOC[idx-1] > soc_target, "wrong direction"
            I_minus_increase = (self.SOC[idx]-soc_target)*self.battery.GetReferenceCapacity()#/self.battery.GetEfficiencyIn()
            #print I_plus_increase
            P_minus_increase = I_minus_increase * self.battery.CalculateDischargeVoltage(self.tau)
            self.B_minus[idx] = 1
            if P_minus_increase > self.battery.GetMaxDelivery(self.tau,self.nu)/self.battery.GetEfficiencyOut() - self.P_out_b[idx]:# "Charge limit exceeded."
                assert idx != 0, "end of the line - no more discharge capacity."
                P_deficit = P_minus_increase - (self.battery.GetMaxDelivery(self.tau,self.nu)/self.battery.GetEfficiencyOut() - self.P_out_b[idx])
                #self.battery.SetSOC(self.SOC[idx-1])
                I_deficit = P_deficit / self.battery.CalculateDischargeVoltage(self.tau)
                new_target = self.SOC[idx-1] - I_deficit/self.battery.GetReferenceCapacity()
                assert new_target <= 1.0+self.eps and new_target >= -1*self.eps, "invalid SOC target. prev soc: "+str(self.SOC[idx-1])+" target:"+str(new_target)
                #print "next match d4:",idx-1,new_target-self.eps
                self.MatchTargetSOC(idx-1,new_target-self.eps)
                #print 'discharge target matched4:',new_target 
                self.ComputeIntervalSOCs()
                return self.GetDischargeDeficit(idx,soc_target)
            return P_minus_increase*self.battery.GetEfficiencyOut(), 0, P_minus_increase
            
    def GetChargeDeficit(self,idx,soc_target):
        """
        based on current SOC and a target, determines the amount of net charge
        to the battery required to match the target SOC.
        idx - period index in interval
        soc_target - target state-of-charge for the battery, e.g., the reset 
            point
        retval - diesel_to_add, delta_P_plus, delta_P_minus
        """
        #print "ch",idx,soc_target, self.SOC[idx]
        if idx == 0:
            self.battery.SetSOC(self.SOC0) 
        else:
            self.battery.SetSOC(self.SOC[idx-1]) 
        if soc_target <= self.SOC[idx]: return 0,0,0
        if self.B_plus[idx] == 1:
            if self.SOC[idx-1] < soc_target: 
                I_plus_increase = (soc_target-self.SOC[idx])*self.battery.GetReferenceCapacity()/self.battery.GetEfficiencyIn()
                P_plus_increase = I_plus_increase * self.battery.CalculateChargeVoltage(self.tau)
                if (P_plus_increase > self.battery.GetMaxCharge(self.tau,self.nu) - self.P_in_b[idx] + self.eps or
                        P_plus_increase > self.RemainingDieselAvailable(idx)+self.eps):# "Charge limit exceeded."
                    assert idx != 0, "end of the line - no more charge capacity. diesel "+str(P_plus_increase > self.RemainingDieselAvailable(idx))+" chargecap"+str(P_plus_increase > self.battery.GetMaxCharge(self.tau,self.nu) - self.P_in_b[idx])
                    P_deficit = max(P_plus_increase - (self.battery.GetMaxCharge(self.tau,self.nu) - self.P_in_b[idx]),
                        P_plus_increase - self.RemainingDieselAvailable(idx)) 
                    #print "P_deficit:",P_deficit,"idx",idx
                    #print "prev SOC:", self.SOC[idx-1]
                    #print "diesel short", P_plus_increase - self.RemainingDieselAvailable(idx)
                    #print "charge short", P_plus_increase - (self.battery.GetMaxCharge(self.tau,self.nu) - self.P_in_b[idx])
                    if idx == 1:
                        self.battery.SetSOC(self.SOC0) 
                    else:
                        self.battery.SetSOC(self.SOC[idx-2])
                    I_deficit = P_deficit / self.battery.CalculateChargeVoltage(self.tau)
                    new_target = self.SOC[idx-1] + self.battery.GetEfficiencyIn()*I_deficit/self.battery.GetReferenceCapacity()
                    assert new_target <= 1.0+self.eps and new_target >= -1*self.eps, "invalid SOC target. prev soc: "+str(self.SOC[idx-1])+" target:"+str(new_target)
                    #print 'new target:',new_target, 'charge1 for period',idx-1
                    self.MatchTargetSOC(idx-1,new_target+self.eps)
                    #print 'charge target matched1:',new_target, 'for period',idx-1 
                    self.ComputeIntervalSOCs()
                    return self.GetChargeDeficit(idx,soc_target)
                if P_plus_increase < 0: 
                    #print "negative P_plus_increase."
                    return 0,0,0
                return P_plus_increase, P_plus_increase, 0
            else:
                I_minus_increase = (self.SOC[idx-1]-soc_target)*self.battery.GetReferenceCapacity()
                P_minus_increase = I_minus_increase * self.battery.CalculateDischargeVoltage(self.tau)
                if P_minus_increase > self.battery.GetMaxDelivery(self.tau,self.nu)/self.battery.GetEfficiencyOut() - self.P_out_b[idx] + self.P_in_b[idx]: #"Discharge limit exceeded."
                    assert idx != 0, "end of the line - no more discharge capacity."
                    P_deficit = P_minus_increase - (self.battery.GetMaxDelivery(self.tau,self.nu)/self.battery.GetEfficiencyOut() - self.P_out_b[idx])
                    self.battery.SetSOC(self.SOC[idx-1])
                    if idx == 1:
                        self.battery.SetSOC(self.SOC0) 
                    else:
                        self.battery.SetSOC(self.SOC[idx-2])
                    I_deficit = P_deficit / self.battery.CalculateDischargeVoltage(self.tau)
                    new_target = self.SOC[idx-1] - I_deficit/self.battery.GetReferenceCapacity()
                    assert new_target <= 1.0+self.eps and new_target >= -1*self.eps, "invalid SOC target. prev soc: "+str(self.SOC[idx-1])+" target:"+str(new_target)
                    #print 'new target:',new_target, 'charge2 for period',idx-1
                    self.MatchTargetSOC(idx-1,new_target-self.eps)
                    #print 'charge target matched2:',new_target 
                    self.ComputeIntervalSOCs()
                    return self.GetDischargeDeficit(idx,soc_target)
                self.B_plus[idx] = 0
                self.B_minus[idx] = 1
                return P_minus_increase*self.battery.GetEfficiencyOut()+self.P_in_b[idx], -1.0*self.P_in_b[idx], P_minus_increase
        elif self.B_minus[idx] == 1:   
            if self.SOC[idx-1] > soc_target: 
                I_minus_reduction = (soc_target-self.SOC[idx])*self.battery.GetReferenceCapacity()
                P_minus_reduction  = I_minus_reduction * self.battery.CalculateDischargeVoltage(self.tau)
                #print I_minus_reduction, soc_target, self.SOC[idx], self.battery.GetReferenceCapacity(), self.battery.CalculateDischargeVoltage(self.tau), self.SOC[idx-1], self.battery.GetSOC()
                #assert P_minus_reduction+self.P_out_b[idx] < self.battery.GetMaxDelivery(self.tau,self.nu)/self.battery.GetEfficiencyOut(), "too much power out in period"+str(idx)
                if P_minus_reduction*self.battery.GetEfficiencyOut() > self.RemainingDieselAvailable(idx):
                    P_deficit = P_minus_reduction*self.battery.GetEfficiencyOut() - self.RemainingDieselAvailable(idx)
                    if idx == 1:
                        self.battery.SetSOC(self.SOC0) 
                    else:
                        self.battery.SetSOC(self.SOC[idx-1])
                    I_deficit = P_deficit / self.battery.CalculateDischargeVoltage(self.tau)
                    new_target = max(self.SOC[idx-1] + I_deficit/self.battery.GetReferenceCapacity(),
                        (I_deficit + self.I_minus[idx-1])*(1+self.battery.GetRateOut())/self.battery.GetReferenceCapacity())
                    assert new_target <= 1.0+self.eps and new_target >= -1*self.eps, "invalid SOC target. prev soc: "+str(self.SOC[idx-1])+" target:"+str(new_target)
                    #print 'new target:',new_target, 'charge3 for period',idx-1
                    self.MatchTargetSOC(idx-1,new_target+self.eps)
                    #print 'charge target matched3:',new_target, 'for period',idx-1
                    self.ComputeIntervalSOCs()
                    return self.GetChargeDeficit(idx,soc_target)
                return P_minus_reduction*self.battery.GetEfficiencyOut(), 0, -1.0*P_minus_reduction
            else:
                I_plus_increase = (soc_target-self.SOC[idx-1])*self.battery.GetReferenceCapacity()/self.battery.GetEfficiencyIn()
                P_plus_increase = I_plus_increase * self.battery.CalculateChargeVoltage(self.tau)
                if (P_plus_increase > self.battery.GetMaxCharge(self.tau,self.nu) - self.P_in_b[idx] or
                        P_plus_increase+self.P_out_b[idx]*self.battery.GetEfficiencyOut() > self.RemainingDieselAvailable(idx)):# "Charge limit exceeded."
                    assert idx != 0, "end of the line - no more charge capacity."
                    P_deficit = max(P_plus_increase - (self.battery.GetMaxCharge(self.tau,self.nu) - self.P_in_b[idx]),
                        P_plus_increase+self.P_out_b[idx]*self.battery.GetEfficiencyOut() - self.RemainingDieselAvailable(idx)) 
                    if idx == 1:
                        self.battery.SetSOC(self.SOC0) 
                    else:
                        self.battery.SetSOC(self.SOC[idx-2])
                    I_deficit = P_deficit / self.battery.CalculateChargeVoltage(self.tau)
                    new_target = self.SOC[idx-1] + self.battery.GetEfficiencyIn()*I_deficit/self.battery.GetReferenceCapacity()
                    assert new_target <= 1.0+self.eps and new_target >= -1*self.eps, "invalid SOC target. prev soc: "+str(self.SOC[idx-1])+" target:"+str(new_target)
                    #print 'new target:',new_target, 'charge4 for period',idx-1
                    self.MatchTargetSOC(idx-1,new_target+self.eps)
                    #print 'charge target matched4:',new_target, 'for period',idx-1 
                    self.ComputeIntervalSOCs()
                    return self.GetChargeDeficit(idx,soc_target)
                if P_plus_increase < 0:
                    #print "negative P_plus_increase."
                    return 0,0,0
                self.B_minus[idx] = 0
                self.B_plus[idx] = 1
                return P_plus_increase+self.P_out_b[idx]*self.battery.GetEfficiencyOut(), P_plus_increase, -1.0*self.P_out_b[idx] 
        else:
            assert self.SOC[idx-1] < soc_target, "wrong direction"
            I_plus_increase = (soc_target-self.SOC[idx])*self.battery.GetReferenceCapacity()/self.battery.GetEfficiencyIn()
            #print I_plus_increase
            P_plus_increase = I_plus_increase * self.battery.CalculateChargeVoltage(self.tau)
            if (P_plus_increase > self.battery.GetMaxCharge(self.tau,self.nu) - self.P_in_b[idx] or
                        P_plus_increase+self.P_out_b[idx]*self.battery.GetEfficiencyOut() > self.RemainingDieselAvailable(idx)):# "Charge limit exceeded."
                    assert idx != 0, "end of the line - no more charge capacity."
                    P_deficit = max(P_plus_increase - (self.battery.GetMaxCharge(self.tau,self.nu) - self.P_in_b[idx]),
                        P_plus_increase+self.P_out_b[idx]*self.battery.GetEfficiencyOut() - self.RemainingDieselAvailable(idx)) 
                    if idx == 1:
                        self.battery.SetSOC(self.SOC0) 
                    else:
                        self.battery.SetSOC(self.SOC[idx-2])
                    I_deficit = P_deficit / self.battery.CalculateChargeVoltage(self.tau)
                    new_target = self.SOC[idx-1] + self.battery.GetEfficiencyIn()*I_deficit/self.battery.GetReferenceCapacity()
                    assert new_target <= 1.0+self.eps and new_target >= -1*self.eps, "invalid SOC target. prev soc: "+str(self.SOC[idx-1])+" target:"+str(new_target)
                    #print 'new target:',new_target, 'charge5 for period',idx-1
                    self.MatchTargetSOC(idx-1,new_target+self.eps)
                    #print 'charge target matched5:',new_target, 'for period',idx-1 
                    self.ComputeIntervalSOCs()
                    return self.GetChargeDeficit(idx,soc_target)
            self.B_plus[idx] = 1
            #assert P_plus_increase < self.battery.GetMaxCharge(self.tau,self.nu) - self.P_in_b[idx], "Charge limit exceeded."
            return P_plus_increase, P_plus_increase, 0
    
    def RemainingDieselAvailable(self,idx):
        return scipy.sum(self.gen_caps - self.P_out_g[idx])
                    
    def AddDieselByPeriod(self,idx,min_soc):
        added_diesel = scipy.sum(self.gen_caps-self.P_out_g,1)
        self.battery.SetSOC(0.0)
        P_plus_est = (-1.0*min_soc*self.battery.GetReferenceCapacity()*
            self.battery.CalculateChargeVoltage(self.inputs["scalars"]["tau"])/
            self.battery.GetEfficiencyIn())
        if P_plus_est > added_diesel[:idx].max():
            per_to_add = added_diesel[:idx].argmax()
            diesel_to_add = added_diesel[:idx].max()
        else:
            per_to_add = self.MinPosDifference(P_plus_est,added_diesel)
            diesel_to_add = P_plus_est
        self.DieselReplaceBattery(diesel_to_add,per_to_add)
            
    def MinPosDifference(self,P_plus_est,added_diesel):
        incumbent=0
        diff=scipy.infty
        surplus = added_diesel-P_plus_est
        for i in range(len(added_diesel)):
            if surplus[i] > 0 and surplus[i] < diff:
                diff = surplus[i]
                incumbent=i
        return incumbent
        
    def AddCharge(self,gen_power,idx):
        """Updates power flow through the battery based on increased power to
        generators.
        gen_power -- increased power to generators"""
        if self.B_minus[idx] == 1:
            if gen_power > self.P_out_b[idx]*self.battery.GetEfficiencyOut():
                self.B_minus[idx] = 0
                self.B_plus[idx] = 1
                self.P_in_b[idx] = gen_power - self.P_out_b[idx]*self.battery.GetEfficiencyOut()
                self.P_out_b[idx] = 0
            else: 
                self.P_out_b[idx] -= gen_power/self.battery.GetEfficiencyOut()
        elif self.B_plus[idx] == 1:
            self.P_in_b[idx] += gen_power
        else:
            self.B_plus[idx] = 1
            self.P_in_b[idx] = gen_power
        
    def DieselReplaceBattery(self,diesel,idx):
        if idx != 0: 
            self.battery.SetSOC(self.SOC[idx-1])
        else:
            self.battery.SetSOC(self.SOC0)
        max_replace_diesel = self.battery.GetMaxCharge(self.tau,self.nu) - self.P_in_b[idx] + self.P_out_b[idx]*self.battery.GetEfficiencyOut()
        diesel = min(max_replace_diesel,diesel)
        #fill existing generators
        for g in range(len(self.G[idx])):
            if self.G[idx][g] == 1:
                surplus = self.gen_caps[g]-self.P_out_g[idx][g]
                if surplus > diesel:
                    self.P_out_g[idx][g] += diesel
                    self.AddCharge(diesel,idx)
                    return None
                else: 
                    self.P_out_g[idx][g] = self.gen_caps[g]
                    self.AddCharge(surplus,idx)
                    diesel -= surplus
        #add generators as needed, smallest that will work first, largest o.w.
        #for g in range(-1,-len(self.gen_caps)-1,-1):
        #    if self.G[idx][g] == 0 and self.gen_caps[g] >= diesel:
        #        self.G[idx][g] = 1
        #        self.P_out_g[idx][g] = diesel
        #        self.AddCharge(diesel,idx)
        #        return None
        ##add the largest generator and try again.
        #for g in range(len(self.gen_caps)):
        #    if self.G[idx][g] == 0:
        #        self.G[idx][g] = 1
        #        self.P_out_g[idx][g] = self.gen_caps[g]
        #        self.AddCharge(self.gen_caps[g],idx)
        #        return self.DieselReplaceBattery(diesel-self.gen_caps[g],idx)    
        
    def SpinningReserveIsMet(self,soc_hat,idx):
        return self.GetSpinningReserve(soc_hat,idx) >= self.reserve[idx]-self.eps
    
    def GetSpinningReserve(self,soc_hat,idx):
        bat_reserve = soc_hat * self.battery.GetEfficiencyOut() * self.battery.GetMaxPower()
        gen_reserve = scipy.sum(self.gen_caps * self.G[idx] - self.P_out_g[idx])
        return bat_reserve + gen_reserve
    
    def AddDieselCapacityOnly(self,deficit,idx):
        assert deficit <= sum(self.gen_caps * (1-self.G[idx])),"IDX"+str(idx)+"\n"+str(self.SOC)+"\nGENCAPS "+str(self.gen_caps)+"\n deficit"+str(deficit)+"\nG "+str(self.G[idx])+"\nP_out_g "+str(self.P_out_g[idx])
        while deficit > self.eps:
            #print 'adddiesel'
            for g in range(-1,-len(self.G[idx])-1,-1):
                if self.G[idx][g] == 0:
                    if self.gen_caps[g] >= deficit or sum(self.G[idx][:g]) == len(self.G[idx][:g]):
                        self.G[idx][g] = 1
                        deficit -= self.gen_caps[g]
                        break
                    
    def CapacityNeeded(self, diesel_to_add, idx):
        return max(0, diesel_to_add - sum(self.gen_caps * self.G[idx]-self.P_out_g[idx]))
    
    def AddDieselOutput(self,diesel_to_add,idx):
        deficit = self.CapacityNeeded(diesel_to_add, idx)
        self.AddDieselCapacityOnly(deficit,idx)
        self.AddGenPower(diesel_to_add,idx)
    
    def ReduceCharge(self,difference):
        self.P_in_b -= difference
        difference = self.ReduceGenPower(difference)
        self.ReducePVOut(difference)
            
    def ReducePVOut(self,difference):
        if difference <= self.P_solar:
            self.P_solar = 0
        else:
            self.P_solar -= difference
            self.reserve -= difference * self.spinning_reserve_factor
    
    def ReduceDieselOutput(self,difference,idx):
        for i in range(-1,-len(self.G[idx])-1,-1):
            if self.G[idx][i] == 1:
                if self.P_out_g[idx][i] <= difference:
                    difference -= self.P_out_g[idx][i] 
                    self.P_out_g[idx][i] = 0.0
                    self.G[idx][i] = 0
                else:
                    self.P_out_g[idx][i] -= difference
                    difference = 0.0
                    break
        return difference
        
    def GetRealSOC(self,idx):
        """gets the actual SOC from power flows given by the IFS, under a 
        nonlinear battery model.  This may be less than zero.
        idx -- value in array that is required
        """
        if idx == 0:
            soc = self.SOC0*1.0
        else:
            soc = self.SOC[idx-1]
        self.battery.SetSOC(soc)
        if self.B_minus[idx] == 1:
            load_delivered = self.P_out_b[idx] * self.battery.GetEfficiencyOut()
            #assert load_delivered > self.battery.GetMaxDelivery(
            #    self.inputs["scalars"]["tau"],self.inputs["scalars"]["nu"]), "P_out exceeds max."
            self.Z_minus[idx], self.I_minus[idx] = self.battery.CalculateRequiredCurrentOut(
                    load_delivered, self.tau)
            soc_hat = (self.battery.GetSOC() - (self.I_minus[idx] * self.tau)/
                    self.battery.GetReferenceCapacity())
            self.Z_plus[idx] = 0
            self.I_plus[idx] = 0
            return soc_hat
        elif self.B_plus[idx] == 1:
            self.Z_plus[idx], self.I_plus[idx] = self.battery.CalculateRequiredCurrentIn(
                    self.P_in_b[idx], self.tau)
            soc_hat = (self.battery.GetSOC() + (self.battery.GetEfficiencyIn() * self.I_plus[idx] * self.tau)/
                    self.battery.GetReferenceCapacity())
            self.Z_minus[idx] = 0
            self.I_minus[idx] = 0
            return soc_hat
        #if B_plus and B_minus are 0, SOC is unchanged.
        return self.battery.GetSOC()
        
        
    #def FillGeneratorPower(self):
    #    """Returns a solution where the generators run at full capacity, 
    #    if able for that time period (based on charge capacity).
    #    """
    #    power_gap = sum(self.G*self.gen_caps-self.P_out_g)
    #    max_charge = self.battery.GetMaxCharge(self.tau,self.nu)
    #    if max_charge < self.P_in_b: 
    #        return None
    #    if power_gap < max_charge + self.P_out_b*self.battery.GetEfficiencyOut() - self.P_in_b:
    #        self.AddGenPower(power_gap)
    #        self.AddBatteryCharge(power_gap)
    #    else:
    #        self.AddGenPower(max_charge + self.P_out_b*self.battery.GetEfficiencyOut() - self.P_in_b)
    #        self.MaxChargeBattery()
            
        
    def AddGenPower(self, added_charge, idx):
        """adds generator power without turning more generators on."""
        reserve = self.G[idx]*self.gen_caps - self.P_out_g[idx]
        for g in range(len(self.G[idx])):
            if reserve[g] > 0.0001:
                if reserve[g] > added_charge:
                    self.P_out_g[idx][g] += added_charge
                    break
                else: 
                    self.P_out_g[idx][g] = self.gen_caps[g]
                    added_charge -= reserve[g]
    
    def MaxChargeBattery(self):
        self.B_minus = 0
        self.B_plus = 1
        self.P_out_b = 0
        self.P_in_b = self.battery.GetMaxCharge(self.tau,self.nu)
    
    def AddBatteryCharge(self,gen_power):
        """Updates power flow through the battery based on increased power to
        generators.
        gen_power -- increased power to generators"""
        if self.B_minus == 1:
            if gen_power > self.P_out_b*self.battery.GetEfficiencyOut():
                self.B_minus = 0
                self.B_plus = 1
                self.P_in_b = gen_power - self.P_out_b*self.battery.GetEfficiencyOut()
                self.P_out_b = 0
            else: 
                self.P_out_b -= gen_power/self.battery.GetEfficiencyOut()
        elif self.B_plus == 1:
            self.P_in_b += gen_power
        else:
            self.B_plus = 1
            self.P_in_b = gen_power
            
    def AddBatteryDischarge(self,power):
        """Updates power flow through the battery based on increased power to
        generators.
        gen_power -- reduced power to generators"""
        if self.B_plus == 1:
            if power > self.P_in_b:
                self.B_plus = 0
                self.B_minus = 1
                self.P_out_b = (power - self.P_in_b)/self.battery.GetEfficiencyOut()
                self.P_in_b = 0
            else: 
                self.P_in_b -= power
        elif self.B_minus == 1:
            self.P_out_b += power/self.battery.GetEfficiencyOut()
        else:
            self.B_minus = 1
            self.P_out_b = power

    def SetSOC(self):
        self.SOC = self.GetRealSOC()*1.0
        assert self.SOC >= 0.0, "negative SOC."
        assert self.SOC <= 1.0, "SOC > 1.0." 
        self.battery.SetSOC(self.SOC*1.0)
        
    
    def StoreDispatch(self,start_period,end_period):
        self.Gs[start_period-1:end_period] = self.G
        self.P_out_gs[start_period-1:end_period] = self.P_out_g
        self.P_out_bs[start_period-1:end_period] = self.P_out_b
        self.P_in_bs[start_period-1:end_period] = self.P_in_b
        self.P_solars[start_period-1:end_period] = self.P_solar
        self.I_pluses[start_period-1:end_period] = self.I_plus
        self.I_minuses[start_period-1:end_period] = self.I_minus
        self.Z_pluses[start_period-1:end_period] = self.Z_plus
        self.Z_minuses[start_period-1:end_period] =self.Z_minus
        self.B_pluses[start_period-1:end_period] = self.B_plus
        self.B_minuses[start_period-1:end_period] = self.B_minus
        self.B_socs[start_period-1:end_period] = self.SOC
        self.F_tildes[start_period-1:end_period] = self.F_tilde
    
    def OutputAllDispatch(self,lead="_new"):
        F_tilde_file = open(self.IFS_dir+"F_tilde"+lead+".csv",'w')
        G_file = open(self.IFS_dir+"G"+lead+".csv",'w')
        P_out_file = open(self.IFS_dir+"P_out"+lead+".csv",'w')
        P_in_file = open(self.IFS_dir+"P_in"+lead+".csv",'w')
        P_solar_file = open(self.IFS_dir+"P_solar"+lead+".csv",'w')
        I_minus_file = open(self.IFS_dir+"I_minus"+lead+".csv",'w')
        I_plus_file = open(self.IFS_dir+"I_plus"+lead+".csv",'w')
        Z_minus_file = open(self.IFS_dir+"Z_minus"+lead+".csv",'w')
        Z_plus_file = open(self.IFS_dir+"Z_plus"+lead+".csv",'w')
        B_minus_file = open(self.IFS_dir+"B_minus"+lead+".csv",'w')
        B_plus_file = open(self.IFS_dir+"B_plus"+lead+".csv",'w')
        B_soc_file = open(self.IFS_dir+"B_soc"+lead+".csv",'w')
        Z_file = open(self.IFS_dir+"Z"+lead+".csv",'w')
        L_file = open(self.IFS_dir+"L"+lead+".csv",'w')
        print "Previous Cost:", self.totalCost
        final_cost = self.CalculateTotalCost()
        print "Final Cost:", final_cost
        print "Increase:", (final_cost-self.totalCost)/self.totalCost
        Z_file.write(self.scenario+"  "+str(final_cost))
        for i in range(len(self.generators)):
            L_file.write(self.scenario+"."+self.generators[i].GetName()+".K"+str(self.generators[i].GetTwin())+"  "+str(scipy.sum(self.Gs,0)[i])+"\n")
        L_file.write(self.scenario+"."+self.battery.GetName()+".K1  "+str(self.GetBatteryLifecycles())+"\n")
        #if t == 1 or t == 5000:
        #    print t, self.G, self.P_out_g
        #    print self.B_minus, self.P_in_b
        #    print self.B_plus
        for t in range(1,self.time_horizon+1):
            for i in range(len(self.generators)):
                if self.Gs[t-1,i] == 1:
                    G_file.write(self.getLead(self.generators[i],t)+"  1\n")
                    P_out_file.write(self.getLead(self.generators[i],t)+"  "+str(self.P_out_gs[t-1,i])+"\n")
            if self.B_pluses[t-1] == 1:
                P_in_file.write(self.getLead(self.battery,t)+"  "+str(self.P_in_bs[t-1])+"\n")
                I_plus_file.write(self.getLead(self.battery,t)+"  "+str(self.I_pluses[t-1])+"\n")
                Z_plus_file.write(self.getLead(self.battery,t)+"  "+str(self.Z_pluses[t-1])+"\n")
                B_plus_file.write(self.getLead(self.battery,t)+"  1\n")
            if self.B_minuses[t-1] == 1:
                P_out_file.write(self.getLead(self.battery,t)+"  "+str(self.P_out_bs[t-1])+"\n")
                I_minus_file.write(self.getLead(self.battery,t)+"  "+str(self.I_minuses[t-1])+"\n")
                Z_minus_file.write(self.getLead(self.battery,t)+"  "+str(self.Z_minuses[t-1])+"\n")
                B_minus_file.write(self.getLead(self.battery,t)+"  1\n")
            B_soc_file.write(self.getLead(self.battery,t)+"  "+str(self.B_socs[t-1])+"\n")
            F_tilde_file.write(self.scenario+".TT"+str(t)+"  "+str(self.F_tildes[t-1])+"\n")
            P_solar_file.write(self.getLead(self.pv_array,t)+"  "+str(self.P_solars[t-1])+"\n")
        G_file.close()
        P_out_file.close()
        P_in_file.close()
        I_minus_file.close()
        I_plus_file.close()
        Z_minus_file.close()
        Z_plus_file.close()
        B_minus_file.close()
        B_plus_file.close()
        Z_file.close()
        L_file.close()
        F_tilde_file.close()
        P_solar_file.close()
        
    def OutputAllDispatchTables(self,lead="_new", output=False):
        F_tilde_file = open(self.IFS_dir+"F_tilde"+lead+".csv",'w')
        G_file = open(self.IFS_dir+"G"+lead+".csv",'w')
        P_out_in_file = open(self.IFS_dir+"P_out_in"+lead+".csv",'w')
        Bat_stats_file = open(self.IFS_dir+"Bat_stats"+lead+".csv",'w')
        Z_file = open(self.IFS_dir+"Z"+lead+".csv",'w')
        L_file = open(self.IFS_dir+"L"+lead+".csv",'w')
        try: final_cost = self.CalculateTotalCost()
        except AttributeError: final_cost = self.totalCost
        if output:
            print "Previous Cost:", self.totalCost
            print "Final Cost:", final_cost
            print "Increase:", (final_cost-self.totalCost)/self.totalCost
        Z_file.write(self.scenario+","+str(final_cost))
        for i in range(len(self.generators)):
            L_file.write(self.scenario+"."+self.generators[i].GetName()+".K"+str(self.generators[i].GetTwin())+","+str(scipy.sum(self.Gs,0)[i])+"\n")
        L_file.write(self.scenario+"."+self.battery.GetName()+".K1,"+str(self.GetBatteryLifecycles())+"\n")
        #if t == 1 or t == 5000:
        #    print t, self.G, self.P_out_g
        #    print self.B_minus, self.P_in_b
        #    print self.B_plus
        P_out_in_file.write("t,")
        for gen in self.generators:
            P_out_in_file.write(gen.GetName()+".K"+str(gen.GetTwin())+",")
        G_file.write("t,")
        for gen in self.generators:
            G_file.write(gen.GetName()+".K"+str(gen.GetTwin())+",")
        G_file.write("\n")
        P_out_in_file.write("Bat_minus_"+self.battery.GetName()+",Bat_plus_"+self.battery.GetName()+",P_Solar\n")
        Bat_stats_file.write("t,I+,Z+,B+,I-,Z-,B-,SOC\n")
        for t in range(1,self.time_horizon+1):
            P_out_in_file.write("TT"+str(t)+",")
            Bat_stats_file.write("TT"+str(t)+",")
            G_file.write("TT"+str(t)+",")
            for i in range(len(self.generators)):
                G_file.write(str(self.Gs[t-1,i])+",")
                P_out_in_file.write(str(self.P_out_gs[t-1,i])+",")
            if self.B_pluses[t-1] == 1:
                P_out_in_file.write("0.0,"+str(self.P_in_bs[t-1])+",")
                Bat_stats_file.write(str(self.I_pluses[t-1])+",")
                Bat_stats_file.write(str(self.Z_pluses[t-1])+",")
                Bat_stats_file.write("1,0,0,0,")
            elif self.B_minuses[t-1] == 1:
                P_out_in_file.write(str(self.P_out_bs[t-1])+",0.0,")
                Bat_stats_file.write("0,0,0,"+str(self.I_minuses[t-1])+",")
                Bat_stats_file.write(str(self.Z_minuses[t-1])+",")
                Bat_stats_file.write("1,")
            else:
                Bat_stats_file.write("0,0,0,0,0,0,")
                P_out_in_file.write("0.0,0.0,")
            Bat_stats_file.write(str(self.B_socs[t-1])+"\n")
            F_tilde_file.write("TT"+str(t)+","+str(self.F_tildes[t-1])+"\n")
            P_out_in_file.write(str(self.P_solars[t-1])+"\n")
            G_file.write("\n")
        G_file.close()
        P_out_in_file.close()
        Bat_stats_file.close()
        Z_file.close()
        L_file.close()
        F_tilde_file.close()
        
    def getLead(self,tech,t):
        if tech.GetType() == "PVArray":
            return self.scenario + "." + tech.GetName() + ".TT" + str(t)
        else:
            return self.scenario + "." + tech.GetName() + ".K" + str(tech.GetTwin()) + ".TT" + str(t)
    
    def RunSolution(self,num_intervals=365,interval_length = 24, output=False):
        for l in range(num_intervals):
            start_period = 1 + l*interval_length
            #print start_period
            end_period = interval_length*(l+1) 
            self.RunDispatchPeriods(start_period,end_period)
        if output:
            try:
                self.OutputAllDispatch(lead="_new")
                self.OutputAllDispatchTables(lead="_nlp",output=True)
            except AttributeError:
                pass

    def RunMIPDispatch(self,t):
        """Runs the MIP Dispatch; used only to calculate cost."""
        self.G = scipy.array([self.GetMIPSolutionValue("G",gen.GetName(),"K"+str(gen.GetTwin()),t) for gen in self.generators])
        self.P_out_g = scipy.array([self.GetMIPSolutionValue("P_out",gen.GetName(),"K"+str(gen.GetTwin()),t) for gen in self.generators])
        try:
            self.P_out_b = self.GetMIPSolutionValue("P_out",self.battery.GetName(),"K1",t)
            self.P_in_b = self.GetMIPSolutionValue("P_in",self.battery.GetName(),"K1",t)
            self.B_plus = self.GetMIPSolutionValue("B_plus",self.battery.GetName(),"K1",t)
            self.B_minus = self.GetMIPSolutionValue("B_minus",self.battery.GetName(),"K1",t)
            self.I_minus = self.GetMIPSolutionValue("I_minus",self.battery.GetName(),"K1",t)
            self.I_plus = self.GetMIPSolutionValue("I_plus",self.battery.GetName(),"K1",t)
            self.Z_minus = self.GetMIPSolutionValue("Z_minus",self.battery.GetName(),"K1",t)
            self.Z_plus = self.GetMIPSolutionValue("Z_plus",self.battery.GetName(),"K1",t)
            self.SOC = self.GetMIPSolutionValue("B_soc",self.battery.GetName(),"K1",t)
        except AttributeError:
            self.P_out_b = 0
            self.P_in_b = 0
            self.B_plus = 0
            self.B_minus = 0
            self.I_minus = 0
            self.I_plus = 0
            self.Z_minus = 0
            self.Z_plus = 0
            self.SOC = 0.0
        self.P_solar = self.GetMIPSolar(t)
        self.F_tilde = self.GetMIPFuel(t)
        self.StoreDispatch(t,t+1)
    

    def RegularSolution(self,output=False):
        for t in range(1,self.time_horizon+1):
            #if t % 1000 == 0: print t
            self.RunMIPDispatch(t)
        if output:
            try:
                self.OutputAllDispatchTables(lead="_old",output=True)
            except AttributeError:
                pass
        
            
    
if __name__ == "__main__":
    local = False    
    import os   
    import sys
    import time
    clock = time.time()
    num_intervals = 365
    #if len(sys.argv) > 1: scenario = sys.argv[1]
    #else: scenario = "ll1"
    outfile = open("MINLP_increaseC5_20p.csv",'w')
    outfile.write("scenario,new_obj,increase,time\n")
    for i in range(1,15):#
        scenario = "ll"+str(i)
        print scenario
        local = True
        if not local: 
            inputs_dir = os.environ['WORK']+"/EEOMC_REPO/OPTIMIZATION/GAMS/"
            IFS_dir = os.environ['WORK']+"/EEOMC_REPO/HEURISTIC/"+scenario+"/"
        else: 
            inputs_dir = "./../OPTIMIZATION/GAMS/"
            IFS_dir = "./OURP/"+scenario+"/"
        M = MINLPSolution(inputs_dir,IFS_dir,scenario)
        #M.SetThreshold(0.00)
        #print "TRYING ORIGINAL SOLUTION."
        #M.RegularSolution(output=True)
        #tc = M.CalculateTotalCost(True)
        print "NEW SOLUTION."
        try:
            M.RunSolution(num_intervals=num_intervals,output=True)
            for i in range(num_intervals):
                if abs(0.5-M.B_socs[(i+1)*24-1]) > 1e-5:
                    print "period", (i+1)*24,"bad SOC", M.B_socs[(i+1)*24-1]
            if M.B_socs.min() < 0:
                print "period", M.B_socs.argmin(), "negative SOC", M.B_socs.min()  
            if M.B_socs.max() > 1.0:
                print "period", M.B_socs.argmax(), "too high SOC", M.B_socs.max()   
            tc = M.CalculateTotalCost(True)
        except AttributeError:
            print "No Battery."
            tc = M.totalCost
        print "TIME:",time.time()-clock
        clock = time.time()
        #tc = M.CalculateTotalCost(True)
        outfile.write(str(i)+","+str(tc)+","+str((tc-M.totalCost)/M.totalCost)+","+str()+"\n")
    