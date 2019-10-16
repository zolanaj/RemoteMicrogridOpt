"""
Alex Zolan (AZ)
FOB Microgrid Optimization Project
National Renewable Energy Laboratory
Created June 9, 2014 (AZ)
Last Updated June 10, 2014 (AZ)

This Class describes the Diesel Generator technology associated with
the FOB Microgrid Optimization project. 

This is part of a heuristic that is expeceted to find an initial feasible
solution to design and dispatch of a system, which will be outputted to a
series of .csv files which will be read by a GAMS or AMPL in order to further
optimize the outcome through a MINLP or MIP solver.
"""

#This defines the class representing the fuel-based Diesel Generator.
#The class defines and gives properties specific to Diesel Generators
#(fuel consumption, lifetime tracking, min/max allowable power ratings,
#min uptime/downtime, and restrictions on ramp up/ramp down of power
#over a single time period).
class DieselGenerator:
    def __init__(self, attributes, twin):#name, twin, purchase_cost, lifecycle_cost, max_life,
                 #efficiency_out, min_power, max_power,
                 #fuel_use_a, fuel_use_b, fuel_use_c, max_ramp_up, max_ramp_down,
                 #min_uptime, min_downtime, tech_max = 6):
        #initialize properties of the generator from the input file.
        self.__attributes = attributes
        self.__generator_type = "DieselGenerator"
        self.__name = attributes["name"]
        self.__twin = twin
        self.__purchase_cost = float(attributes["c_tilde"])
        self.__lifecycle_cost = float(attributes["epsilon"])
        self.__max_life = float(attributes["lj_max"])
        self.__efficiency_out = float(attributes["eta_out"])
        self.__min_power = float(attributes["p_min"])
        self.__max_power = float(attributes["p_max"])
        if "ag" in attributes.keys(): self.__fuel_use_a = float(attributes["ag"])
        else: self.__fuel_use_a = 0.0
        self.__fuel_use_b = float(attributes["bg"])
        self.__fuel_use_c = float(attributes["cg"])
        if "r_up" in attributes.keys(): self.__max_ramp_up = float(attributes["r_up"])
        else: self.__max_ramp_up = 1.0
        if "r_down" in attributes.keys(): self.__max_ramp_down = float(attributes["r_down"])
        else: self.__max_ramp_down = 1.0
        if "t_up" in attributes.keys(): self.__min_uptime = float(attributes["t_up"])
        else: self.__min_uptime = 1.0
        if "t_down" in attributes.keys(): self.__min_downtime = float(attributes["t_down"])
        else: self.__min_downtime = 1.0
        if "tech_max" in attributes.keys(): self.__tech_max = float(attributes["tech_max"])
        else: self.__tech_max = 3  #default.  This should be reset.
        if "w_tilde" in attributes.keys(): self.__weight = float(attributes["w_tilde"])
        else: self.__weight = 0.0
        if "volume" in attributes.keys(): self.__weight = float(attributes["volume"])
        else: self.__volume = 0.0
        #we need to track cycles used, initialize at max lifetime.
        self.__life_remaining = self.__max_life
        #assume any generator is idle and ready to start up at the time of
        #initialization.  Also, assume it is not selected for dispatch yet.
        self.__is_running = False
        self.__current_uptime = 0
        self.__current_downtime = self.__min_downtime
        self.__is_dispatched = True
        #we need a reference of the most recent power output, for ramp-up
        #and ramp-down constraints.
        self.__last_power_out = 0.0
        #self.__fuel_cost = 0
        #track how much fuel has been used so far.
        self.__fuel_consumed = 0.0
        self.__fuel_cost = 0.0
        self.__variable_cost = 0.0
        self.__fuel_curve = []
        
        
    #define accessors for some of the key properties.
    def GetType(self): return self.__generator_type
    def GetName(self): return self.__name
    def GetTwin(self): return self.__twin
    def GetPurchaseCost(self): return self.__purchase_cost
    def GetLifeCycleCost(self): return self.__lifecycle_cost
    def GetMaxLife(self): return self.__max_life
    def GetLifeRemaining(self): return self.__life_remaining
    def GetEfficiencyOut(self): return self.__efficiency_out
    def GetFuelUseA(self): return self.__fuel_use_a
    def GetFuelUseB(self): return self.__fuel_use_b
    def GetFuelUseC(self): return self.__fuel_use_c
    def GetMinPower(self): return self.__min_power
    def GetMaxPower(self): return self.__max_power
    def GetMinUptime(self): return self.__min_uptime
    def GetMinDowntime(self): return self.__min_downtime
    def GetMaxRampUp(self): return self.max_ramp_up
    def GetMaxRampDown(self): return self.max_ramp_down
    def GetCurrentUptime(self): return self.__current_uptime
    def GetCurrentDowntime(self): return self.__current_downtime
    def GetVariableCost(self): return self.__variable_cost
    def IsDispatched(self): return self.__is_dispatched
    def GetTechMax(self):  return self.__tech_max
    def GetFuelCost(self): return self.__fuel_cost
    def GetFuelUsed(self): return self.__fuel_consumed
    def GetWeight(self): return self.__weight
    def GetVolume(self): return self.__volume

    #define a mutator for fuel_cost.
    #def SetFuelCost(self, fuel_cost): self.__fuel_cost = fuel_cost
    def SetTechMax(self, tech_max): self.__tech_max = tech_max
    def GetFuelCurve(self): return self.__fuel_curve

    #define mutators that select or deselect generators for dispatch.
    def SelectForDispatch(self): self.__is_dispatched = True
    def DeselectForDispatch(self): self.__is_dispatched = False
    
    def ReduceFuelCost(self,factor):
        for i in range(len(self.__fuel_cost)):
            self.__fuel_cost /= factor

    #This determines whether the generator can be run in the next time period.
    def IsAvailable(self,steplength, nu): return (self.__is_running or
                                   self.GetCurrentDowntime() >=
                                   self.GetMinDowntime()) and \
                                   self.__life_remaining > steplength * nu
    
    #This determines whether the generator can be off in the current time period.
    def CanBeStopped(self): return not (self.__is_running and self.GetCurrentUptime <
                                        self.GetMinUptime)
                                    

    #This function returns the amount of power than can be delivered to meet
    #load, net of the efficiency (eta_in, not fuel use) of the generator.
    def GetPowerAvailable(self): return self.__efficiency_out * self.__max_power

    """
    #This calculates the fuel usage over a single time period giiven the
    #fuel cooefficients, length off time, and power generated.
    def AddFuelUsed(self, power_out, steplength):
        fuel_consumed = 
        self.__fuel_consmed += fuel_consumed
        return fuel_consumed
    """        
    
    #This function determines the amount of fuel that is burned based on 
    #the power output of the generator.  It uses a fuel curve to find 
    #this.
    def GetFuelConsumed(self, power_out, steplength):
        return steplength * (self.__fuel_use_a*pow(power_out,2)/1000000.0
                         + self.__fuel_use_b * power_out/1000.0 
                         + self.__fuel_use_c * 1.0 ) 
    
    #This function runs the generator for a single time period and adds the fuel
    #consumed over the time period steplength.  Idling is simply keeping the
    #generator running at zero watts out and is covered in this function.
    def RunDieselGenerator(self, power_out, fuel_cost, steplength, nu):
        assert self.IsAvailable(steplength,nu), \
            "Generator %s cannot be run - insufficient downtime." % self.__name
        if self.__is_running:
            self.__current_uptime += steplength
        else:
            self.__is_running = True
            self.__current_uptime = steplength
            self.__current_downtime = 0
        fuel_consumed = self.GetFuelConsumed(power_out, steplength)
        self.__fuel_consumed += fuel_consumed
        self.__variable_cost += self.__lifecycle_cost * steplength
        self.__fuel_cost += fuel_consumed * fuel_cost
        self.__life_remaining -= steplength * nu
        return fuel_consumed

    #This function stops the running of the Diesel Generator, by editing the
    #is_running state (if needed) and adjusting current downtime.
    def StopDieselGenerator(self, steplength):
        assert (not self.__is_running or self.__current_uptime >= self.__min_uptime),\
            "Generator %s cannot be stopped - insufficient uptime." % self.__name
        if not self.__is_running:
            self.__current_downtime += steplength
        else:
            self.__is_running = False
            self.__current_downtime = steplength
            self.__current_uptime = 0
            
    """
    #Return the total variable cost of dispatch.  This is the cost of
    #Generator depreciation and fuel usage.
    def CalculateVariableCost(self):
        return self.__fuel_cost * self.__fuel_consumed + self.__lifecycle_cost*(
            self.__max_life - self.__life_remaining)
    """
    
    #This function will create a copy of the Diesel Generator, setting the
    #twin number as the number given.
    def CopyDieselGenerator(self,twin):
        gen = DieselGenerator(
            self.__attributes,
            twin
            )
        #gen.SetFuelCost(self.__fuel_cost)
        return gen

