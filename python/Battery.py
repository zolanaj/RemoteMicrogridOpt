"""
Alex Zolan (AZ)
FOB Microgrid Optimization Project
National Renewable Energy Laboratory
Created June 9, 2014 (AZ)
Last Updated June 10, 2014 (AZ)

This Class describes the Energy Storage (Battery) technology associated
with the FOB Microgrid Optimization project.

This is part of a heuristic that is expected to find an initial feasible
solution to design and dispatch of a system, which will be outputted to a
series of .csv files which will be read by a GAMS or AMPL in order to further
optimize the outcome through a MINLP or MIP solver.

Update 8/14/2014: Battery now has maximum and minimum state of charge that can
serve as upper/lower bounds in the optimization program.
  
Update 9/10/2014: Battery lifetime now in terms of full cycles, rather than 
change of direction.  Lifetime now reduced as charges and discharges are made.

Update 9/10/2014: Minimum Current added.  
"""

#This defines the Battery Class.  Its main components are voltage rating,
#reference capacity and related capacity factors, max power rating,
#internal resistance and duty cycle constraints.  Dispatch instructions are
#also included here.
class Battery:
    def __init__(self, attributes, twin):#name, twin, purchase_cost, lifecycle_cost, max_life, 
                #efficiency_in, efficiency_out, min_power, max_power, resistance, voltage_a,
                #voltage_b, reference_capacity, max_rate_out, max_rate_in,
                #duty_cycle_min, duty_cycle_max, duty_voltage, a_soc, d_soc,
                #a_life_temp, b_life_temp, c_life_temp,tech_max,min_current = 0):
        #initialize properties of the battery from the input file.
        self.__generator_type = "Battery"
        self.__name = attributes["name"]
        self.__twin = twin
        self.__purchase_cost = float(attributes["c_tilde"])
        self.__lifecycle_cost = float(attributes["epsilon"])
        self.__max_life = float(attributes["lj_max"])
        self.__efficiency_in = float(attributes["eta_in"])
        self.__efficiency_out = float(attributes["eta_out"])
        self.__min_power = float(attributes["p_min"])
        self.__max_power = float(attributes["p_max"])
        self.__resistance = float(attributes["r_int"])
        self.__voltage_a = float(attributes["a_v"])
        self.__voltage_b = float(attributes["b_v"])
        self.__reference_capacity = float(attributes["c_ref"])
        self.__c_out = float(attributes["c_out"])
        self.__c_in = float(attributes["c_in"])
        self.__avg_current = self.__reference_capacity  #PLACEHOLDER - ADD TO INPUT FILE
        self.__a_soc = float(attributes["a_soc"])
        self.__d_soc = float(attributes["d_soc"])
        if "weight" in attributes.keys(): self.__weight = float(attributes["weight"])
        else: self.__weight = 0.0
        if "volume" in attributes.keys(): self.__weight = float(attributes["volume"])
        else: self.__volume = 0.0
        if "i_min" in attributes.keys(): self.__min_current = float(attributes["i_min"])
        else: self.__min_current = 0
        if "a_l" in attributes.keys(): self.__a_life_temp = float(attributes["a_l"])
        if "b_l" in attributes.keys(): self.__b_life_temp = float(attributes["b_l"])
        if "c_l" in attributes.keys(): self.__c_life_temp = float(attributes["c_l"])
        if "tech_max" in attributes.keys(): self.__tech_max = int(attributes["tech_max"])
        else: self.__tech_max = 1
        #initialize the life remaining (counting for variable cost calculation)
        self.__life_remaining = self.__max_life * 1.0
        #initialize the state of charge, its bounds, and the variable cost incurred.
        self.__soc = 0.0
        self.__min_soc = 0.0
        self.__max_soc = 1.0
        self.__variable_cost = 0.0
        
        #indicator for whether the last active dispatch on the battery was a
        #charge - used for counting cycles.
        self.__last_charged = True
        #indicator for the state of charge as of the most recent charge (for
        #depth of discharge calculation)
        self.__last_charge_soc = 0.0
        
    #define accessors for some of the key properties.
    def GetType(self): return self.__generator_type
    def GetName(self): return self.__name
    def GetTwin(self): return self.__twin
    def GetPurchaseCost(self): return self.__purchase_cost
    def GetLifeRemaining(self): return self.__life_remaining
    def GetLifeCycleCost(self): return self.__lifecycle_cost
    def GetMaxLife(self): return self.__max_life
    def GetSOC(self): return self.__soc
    def GetResistance(self): return self.__resistance
    def GetVoltageA(self): return self.__voltage_a
    def GetVoltageB(self): return self.__voltage_b
    def GetAvgCurrent(self): return self.__avg_current
    def GetEfficiencyIn(self): return self.__efficiency_in
    def GetEfficiencyOut(self): return self.__efficiency_out
    def GetReferenceCapacity(self): return self.__reference_capacity
    def GetMinCurrent(self): return self.__min_current
    def GetVariableCost(self): return self.__variable_cost
    def GetTechMax(self):  return self.__tech_max
    def GetASOC(self):  return self.__a_soc
    def GetDSOC(self):  return self.__d_soc
    def GetMaxPower(self):  return self.__max_power
    def GetMinPower(self):  return self.__min_power
    def GetRateOut(self): return self.__c_out
    def GetRateIn(self): return self.__c_in
    def GetIUMinus(self): return self.__reference_capacity / (1+self.__c_out)
    def GetIUPlus(self): return self.__reference_capacity / self.__c_in
    def GetMinSOC(self): return self.__min_soc
    def GetMaxSOC(self): return self.__max_soc
    
    #define mutators
    def SetSOC(self, soc): self.__soc = soc
    def SetMinSOC(self, soc): self.__min_soc = soc
    def SetMaxSOC(self, soc): self.__max_soc = soc
    def SetLastCharged(self, last_charged): self.__last_charged = last_charged
    def SetLastChargeSOC(self, soc): self.__last_charge_soc = soc
    def SetMinCurrent(self, min_current): self.__min_current = min_current
    def SetTechMax(self, tech_max): self.__tech_max = tech_max
    
    #This function calculates the voltage if the battery is to be discharged.
    def CalculateDischargeVoltage(self, steplength):
        return self.__voltage_a * self.__soc + self.__voltage_b - (
            self.__avg_current * self.__resistance ) 
    
    #This function calculates the voltage if the battery is to be charged.
    def CalculateChargeVoltage(self, steplength):
        return self.__voltage_a * self.__soc + self.__voltage_b + (
            self.__avg_current * self.__resistance )
    
    def CalculateIdleVoltage(self):
        return self.__voltage_a * self.__soc + self.__voltage_b
    
    #This function calculates the effective capacity of the battery.
    def CalculateAdjustedCapacity(self, current_out):
        return self.__reference_capacity - (self.__c_out *
                                            current_out)
    
    #This function returns the minimum power the battery can deliver in a given
    #time period, net of efficiency.
    def GetMinDelivery(self,steplength):
        return max(self.__min_power,(self.CalculateDischargeVoltage(steplength)*self.__min_current))* self.__efficiency_out
                

    def GetMinCharge(self,steplength):
        return max(((self.CalculateChargeVoltage(steplength)*self.__min_current)
                / self.__efficiency_in),self.__min_power)

    #This function returns the maximum power the battery can deliver in a given
    #time period, net of efficiency.
    def GetMaxDelivery(self, steplength, nu): 
        power_delivered = min(min(
            (self.__reference_capacity * self.__soc) /
            (steplength + self.__c_out),
            (self.__reference_capacity * (self.__soc - self.__min_soc)),
            ) * self.CalculateDischargeVoltage(steplength),self.__max_power) * self.__efficiency_out
        #if self.__soc < self.__min_soc or self.__soc > self.__max_soc: print self.__soc
        if self.CanBeDischarged(power_delivered, steplength, nu):
            return power_delivered
        else: return 0
    """    
    def CalculateSpinningReserve(self, load):  
        if load < 0:
            return (
        self.__efficiency_out * self.__max_power * self.__soc 
        )
    """    
    #This function returns the maximum power (W) that can be sent to charge
    #the battery.
    def GetMaxCharge(self, steplength, nu):
        max_power_in = min(self.__max_power,
            self.CalculateChargeVoltage(steplength) * 
             min((self.__max_soc - self.__soc) *
             self.__reference_capacity / (self.__efficiency_in*steplength),
             steplength * self.__reference_capacity / (self.__c_in)
            ))  
        if self.CanBeCharged(max_power_in, steplength, nu):
            return max_power_in
        else: return 0
    
    #This function determines the current necessary to deliver a given amount of
    #power to be delivered to serve load.
    def CalculateRequiredCurrentOut(self,power_delivered,steplength):
        I_minus = power_delivered/(self.__efficiency_out*
                self.CalculateDischargeVoltage(steplength)  )
        Z_minus = I_minus * self.__soc
        return Z_minus, I_minus

    #This function determines the current necessary to deliver a given amount of
    #power to be sent to the battery.
    def CalculateRequiredCurrentIn(self,power_delivered,steplength):
        I_plus = power_delivered / (
                self.CalculateChargeVoltage(steplength) )
        Z_plus = I_plus * self.__soc
        return Z_plus, I_plus

    #This function assesses whether or not the battery can be discharged at the
    #rate of power delivered.
    def CanBeDischarged(self, power_delivered, steplength, nu):
        #Z_minus, I_minus = self.CalculateRequiredCurrentOut(power_delivered,
        #                                               steplength)
        #return self.hasLifeRemaining(I_minus, steplength, nu)
        return True

    #This function assesses whether or not the battery can be charged at the
    #rate of power delivered.
    def CanBeCharged(self, power_delivered, steplength, nu):
        Z_plus, I_plus = self.CalculateRequiredCurrentIn(power_delivered,
                                                 steplength)
        #capacity = self.__reference_capacity 
        return self.hasLifeRemaining(Z_plus, I_plus,steplength, nu) 
        #(cur_in  <= capacity * steplength / self.__c_in
        #         and (cur_in <= capacity * (1-self.__soc) 
        #         / steplength ) and 
    
    #This function denotes the transition from charging to discharging the
    #battery. Under the current definition, we track only change of direction
    #and the state of charge at the last charge.  
    #(No longer needed as of 8/14/2014)
    def DischargeCycle(self):
        self.__last_charged = False
        self.__last_charged_soc = self.__soc
        self.__life_remaining -= 1
    
    #This function denotes a transition from discharging to charging the
    #battery.  Under the current definition of lifecycle tracking
    #we are only looking at changes in direction of current, so we only
    #change the related boolean here.
    def ChargeCycle(self): self.__last_charged = True
    
    #Update the life left in the battery and record the     
    def UpdateLifetime(self, Z_plus, I_plus, steplength, nu):
        life_expended = steplength * (I_plus * self.__a_soc - 
                (Z_plus * self.__d_soc) ) / ( 
                        self.__reference_capacity)
        self.__life_remaining -= life_expended * nu
        self.__variable_cost += self.__lifecycle_cost * life_expended
        return life_expended
        
    #This function discharges the battery, asserting that first it is feasible
    #to do so, and then by updating state of charge and lifecycles if needed.
    #Return P_out and L_bkt (lifetime lost)
    def DischargeToLoad(self, power_delivered, steplength, nu):
        #add an assertion on the adjusted capacity being feasible given the
        #power to be delivered
        #assert self.CanBeDischarged(power_delivered, steplength), \
        #       "Cannot Discharge Battery %r" % self.__name
        #if self.__last_charged:
        #    self.DischargeCycle()
        Z_minus, I_minus = self.CalculateRequiredCurrentOut(power_delivered,
                                                       steplength)
        #print "Starting SOC:", self.__soc
        #print "Outgoing current:", current_out
        #life_expended = self.UpdateLifetime(current_out, steplength, nu) #Lifetime is updated with new SoC 
        self.SetSOC(self.__soc - (I_minus * steplength)/
                    self.__reference_capacity)
        #print "Final SOC:", self.__soc
        #assert  0 <= self.__soc <= 1, "%r soc out of range, %f" % self.__name % self.__soc
        return Z_minus, I_minus, power_delivered/self.__efficiency_out
    
    def hasLifeRemaining(self, Z_plus, I_plus, steplength, nu):
        """asserts that the battery is able to be used."""
        return self.__life_remaining >= nu * (steplength/
                (2*self.__reference_capacity)) * \
                (self.__a_soc*I_plus-self.__d_soc*Z_plus) 
    
    #This function charges the battery, asserting that first it is feasible
    #to do so, and then by updating state of charge and lifecycles if needed.
    #Return P_in.
    def ChargeToBattery(self, power_delivered, steplength, nu):
        #add an assertion on the adjusted capacity being feasible given the
        #power to be delivered
        #assert self.CanBeCharged(power_delivered, steplength), \
        #       "Cannot Charge Battery %r" % self.__name
        #if not self.__last_charged:
        #    self.ChargeCycle()
        Z_plus, I_plus  = self.CalculateRequiredCurrentIn(power_delivered,
                                                     steplength)
        """        
        if self.GetTwin() == 1:                                             
            print "Starting SOC:", self.__soc
            print "Power Delivered:", power_delivered
            print "Voltage:", self.CalculateChargeVoltage(steplength)
            print "Incoming current:", current_in  
        """
        life_expended = self.UpdateLifetime(Z_plus, I_plus, steplength, nu)                                        
        self.SetSOC(self.__soc + self.__efficiency_in * 
                    (I_plus * steplength) / self.__reference_capacity)  
        #Lifetime is updated with new SoC
        #if self.GetTwin() == 1: print "Final SOC:", self.__soc
        #assert  0 <= self.GetSOC() <= 1, "soc above range, %r" % self.__name % str(self.__soc)     
        return Z_plus, I_plus, power_delivered, life_expended

    #This function creates a copy of the Battery, with the twin number set as
    #given.
    def CopyBattery(self,twin):
        bat = Battery(
            {"name":self.__name,
            "c_tilde":self.__purchase_cost,
            "epsilon":self.__lifecycle_cost,
            "lj_max":self.__max_life,
            "eta_in":self.__efficiency_in,
            "eta_out":self.__efficiency_out,
            "p_min":self.__min_power,
            "p_max":self.__max_power,
            "r_int":self.__resistance,
            "a_v":self.__voltage_a,
            "b_v":self.__voltage_b,
            "c_ref":self.__reference_capacity,
            "c_out":self.__c_out,
            "c_in":self.__c_in,
            #"d_min":self.__duty_cycle_min,
            #"d_max":self.__duty_cycle_max,
            #"v_conv":self.__duty_voltage,
            "a_soc":self.__a_soc,
            "d_soc":self.__d_soc,
            #"a_l":self.__a_life_temp, 
            #"b_l":self.__b_life_temp, 
            #"c_l":self.__c_life_temp,
            "tech_max":self.__tech_max,
            "i_min":self.__min_current
            },
            twin
            )
        bat.SetSOC(self.__soc)
        bat.SetLastCharged(self.__last_charged)
        bat.SetLastChargeSOC(self.__last_charge_soc)
        bat.SetMinSOC(self.__min_soc)
        bat.SetMaxSOC(self.__max_soc)
        return bat
