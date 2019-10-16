"""
This class describes the Photovolitaic Array (PVArray) technology associated
with the FOB Microgrid Optimization project. This is part of a heuristic that 
is expected to find an initial feasible
solution to design and dispatch of a system, which will be outputted to a
series of .csv files which will be read by a GAMS or AMPL in order to
optimize the outcome through a MINLP or MIP solver.
"""

#This defines the Photovoltaic (PV) array.  Its main components are efficiency
#and surface area, and it will return power output available based on these
#factors.  Because the power available (in Watts per meter squared) is going
#to be available as a standalone parameter that is net of everything else, there
#is no need to include all the PVWatts parameters in the scope of this diispatch
#algorithm.
import scipy
class PVArray:
    def __init__(self, attributes, twin):
        #initialize properties of the PV array from the input file.
        #note that power rating for solar is in kW.
        self.__attributes = attributes
        self.__generator_type = "PVArray"
        self.__name = attributes["name"]
        self.__twin = twin
        self.__purchase_cost = attributes["purchase_cost"]
        #self.__efficiency_out = attributes["eta_out"]
        self.__surface_area = attributes["footprint"]
        self.__power_rating = attributes["powerRating"]
        self.__pv_availability = {}
        self.__tech_max = 75
        self.__num_panels = 1
        if "weight" in attributes.keys(): self.__weight = float(attributes["weight"])
        else: self.__weight = 0.0
        if "volume" in attributes.keys(): self.__volume = float(attributes["volume"])
        else: self.__volume = 0.0
        
    #define accessors for some of the key properties.
    def GetType(self): return self.__generator_type
    def GetName(self): return self.__name
    def GetTwin(self): return self.__twin
    def GetPurchaseCost(self): return self.__purchase_cost
    def GetEfficiencyOut(self): return self.__efficiency_out
    def GetPowerRating(self): return self.__power_rating
    def GetPVAvailability(self): return self.__pv_availability
    def GetSurfaceArea(self): return self.__surface_area
    def GetTechMax(self):  return self.__tech_max
    def GetNumberOfPanels(self): return self.__num_panels
    #Return the total variable cost of dispatch.  This is zero for a PV array.
    def GetVariableCost(self): return 0.0
    def GetWeight(self): return self.__weight
    def GetVolume(self): return self.__volume    
        
    #define mutator for PV availability - this is used by the PVWattsRun module
    def SetPVAvailability(self, pv_availability):
        #print pv_availability.keys() 
        if type(pv_availability) == list:
            self.__pv_availability = scipy.array(pv_availability)
        else: self.__pv_availability = pv_availability
                
    def ConvertPVAvailability(self, time_horizon):
        """converts dict to array"""
        arr = scipy.zeros(time_horizon)
        for key in self.__pv_availability.keys():
            t = int(key[2:])
            arr[t-1] = self.__pv_availability[key] 
        self.__pv_availability = arr
        #print self.__pv_availability[0:24]      
    
    def SetTechMax(self, tech_max): self.__tech_max = tech_max
        
    def SetNumberOfPanels(self, num_panels):
        self.__num_panels = num_panels
        self.SetPowerRating(self.__power_rating * num_panels)
        
    #define mutator for surface area: used now that the solution is a
    #continuous variable
    def SetPowerRating(self,power_rating): 
        ratio = power_rating / self.__power_rating
        self.__surface_area *= ratio
        self.__power_rating *= ratio
        self.__purchase_cost *= ratio
        self.__weight *= ratio
        self.__volume *= ratio
        self.__pv_availability = ratio * self.__pv_availability
    
    #This calculates the power available through PV, given weather data.
    def CalculatePowerAvailable(self, t): 
        return self.__pv_availability[t] 

    #Makes a copy of the PVArray, with the twin and size number as given.
    def CopyPVArray(self,twin,num_panels):
        array = PVArray(
            self.__attributes,
            twin
            )
        array.SetPVAvailability(self.__pv_availability)
        array.SetNumberOfPanels(num_panels)
        return array
