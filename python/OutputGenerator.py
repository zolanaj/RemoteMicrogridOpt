"""
Alex Zolan (AZ)
FOB Microgrid Optimization Project
National Renewable Energy Laboratory
Last Updated May 4, 2015 (AZ)

This Class generates solutions based on the results of the algorithm.
The heuristic gives us a lot of continuous variables - P_out, P_in,
F_tilde, etc.  However, an initial solution in GAMS is required to
have all the solutions present.

This is part of a heuristic that is expeceted to find an initial feasible
solution to design and dispatch of a system, which will be outputted to a
series of .csv files which will be read by a GAMS or AMPL in order to further
optimize the outcome through a MINLP or MIP solver.
"""

class OutputGenerator:
    def __init__(self,dispatcher,results,soc_start):
        self.__dispatcher = dispatcher
        self.__results = results
        self.__threshold = 0.01 #change to detect for binaries.
        self.__soc_start = soc_start

    def CalculateBatteryBinaries(self):
        B_soc = self.__results['B_soc']
        P_out_b = self.__results['P_out_b']
        P_in_b = self.__results['P_in_b']
        #Initiate variable arrays
        B_plus = [[1 if P_in_b[t][b] > self.__threshold else 0 for b in range(len(self.__dispatcher.GetBatteries()))] for t in range(len(B_soc))]
        B_minus = [[1 if P_out_b[t][b] > self.__threshold else 0  for b in range(len(self.__dispatcher.GetBatteries()))] for t in range(len(B_soc))]
        #Z_plus = self.__results['Z_plus']
        #Z_minus = self.__results['Z_minus']
        #I_plus = self.__results['I_plus']
        #I_minus = self.__results['I_minus']
        #Y_plus = I_plus
        #Y_minus = I_minus
        return {"B_plus":B_plus, "B_minus":B_minus}

    def CalculateVoltCurrentCap(self, battery_binaries):
        soc_start = self.__soc_start
        a_v = [bat.GetVoltageA() for bat in self.__dispatcher.GetBatteries()]
        b_v = [bat.GetVoltageB() for bat in self.__dispatcher.GetBatteries()]
        r_int = [bat.GetResistance()
                     for bat in self.__dispatcher.GetBatteries()]
        i_avg = [bat.GetAvgCurrent()
                     for bat in self.__dispatcher.GetBatteries()]
        #c_ref = [bat.GetReferenceCapacity()
        #             for bat in self.__dispatcher.GetBatteries()]
        B_soc = self.__results["B_soc"]
        B_plus = battery_binaries["B_plus"]
        B_minus = battery_binaries["B_minus"]
        times = range(len(B_soc))
        bats = range(len(B_soc[0]))
        #generate all steps from t=2 on (including dummy values for t=1)
        V_soc = [[a_v[bat]*B_soc[t-1][bat] + b_v[bat] + 
                    r_int[bat]*i_avg[bat]*(B_plus[t][bat]-B_minus[t][bat]) 
                    if ( (B_plus[t][bat]+B_minus[t][bat]) >= 1) 
                    else a_v[bat]*B_soc[t-1][bat] for bat in bats] for t in times]
        #now generate the first step voltage, using soc_start
        V_soc[0] = [a_v[bat]*soc_start + b_v[bat] + 
                    r_int[bat]*i_avg[bat]*(B_plus[0][bat]-B_minus[0][bat]) 
                    if ( (B_plus[0][bat]+B_minus[0][bat]) >= 1) 
                    else a_v[bat]*B_soc[0][bat] for bat in bats]
        #generate currents using P_in/P_out and V_soc
        #I_minus = [[P_out[t][bat]/V_soc[t][bat] if V_soc[t][bat] != 0 else 0 for bat in bats] for t in times]
        #I_plus = [[P_in[t][bat]/V_soc[t][bat] if V_soc[t][bat] != 0 else 0 for bat in bats] for t in times]
        return V_soc#, I_minus, I_plus

    def CalculateGeneratorBinaries(self):
        G = self.__results["G"]
        gens = range(len(G[0]))
        times = range(len(G))
        G_start = [[1 if G[t][gen] == 1 and G[t-1][gen] == 0 else 0 for gen in gens] for t in times]
        G_start[0] = G[0]
        G_stop = [[1 if G[t][gen] == 0 and G[t-1][gen] == 1 else 0 for gen in gens] for t in times]
        G_stop[0] = [0 for gen in gens]
        return G_start, G_stop 

    def GenerateOutputFiles(self,scenario_name = None,output_lead = None,append=False):
        if append: filewrite = 'a'
        else: filewrite = 'w'
        #generate additional battery and generator parameters
        battery_binaries = self.CalculateBatteryBinaries()
        V_soc = self.CalculateVoltCurrentCap(battery_binaries)
        G_start, G_stop = self.CalculateGeneratorBinaries()
        gen_names = [gen.GetName()+'.K'+str(gen.GetTwin()) for gen in self.__dispatcher.GetDieselGenerators()]
        pv_names = [pv.GetName() for pv in self.__dispatcher.GetPVArrays()]
        pv_name_panels = [pv.GetName()+'   '+str(pv.GetNumberOfPanels()) for pv in self.__dispatcher.GetPVArrays()]
        bat_names = [bat.GetName()+'.K'+str(bat.GetTwin()) for bat in self.__dispatcher.GetBatteries()]
        #Set up variables
        B_soc = self.__results["B_soc"]
        P_out_b = self.__results["P_out_b"]
        P_in_b = self.__results["P_in_b"]
        I_minus = self.__results["I_minus"]
        I_plus = self.__results["I_plus"]
        P_out_s = self.__results["P_out_s"]
        P_out_g = self.__results["P_out_g"]
        G = self.__results["G"]
        L_gen = [sum([G[elt][gg] for elt in range(len(G))]) for gg in range(len(G[0]))]
        L = self.__results["L"]
        #L_bkt = self.__results["L_bkt"]
        F_tilde = self.__results["F_tilde"]
        B_plus = battery_binaries["B_plus"]
        B_minus = battery_binaries["B_minus"]
        #Y_plus = I_plus
        #Y_minus = I_minus
        Z_plus = self.__results["Z_plus"]
        Z_minus = self.__results["Z_minus"]     
        #Set up output files
        P_out_file = open(output_lead+'P_out.csv',filewrite) 
        #P_out_file.write(',')
        P_solar_file = open(output_lead+'P_solar.csv',filewrite) 
        P_in_file = open(output_lead+'P_in.csv',filewrite)
        F_tilde_file = open(output_lead+'F_tilde.csv',filewrite)
        B_soc_file = open(output_lead+'B_soc.csv',filewrite)
        B_plus_file = open(output_lead+'B_plus.csv',filewrite)
        B_minus_file = open(output_lead+'B_minus.csv',filewrite)
        G_file = open(output_lead+'G.csv',filewrite)
        G_start_file = open(output_lead+'G_start.csv',filewrite)
        G_stop_file = open(output_lead+'G_stop.csv',filewrite)
        V_soc_file = open(output_lead+'V_soc.csv',filewrite)
        I_plus_file = open(output_lead+'I_plus.csv',filewrite)
        I_minus_file = open(output_lead+'I_minus.csv',filewrite)
        Z_plus_file = open(output_lead+'Z_plus.csv',filewrite)
        Z_minus_file = open(output_lead+'Z_minus.csv',filewrite)
        Z_file = open(output_lead+'Z.csv',filewrite)
        Z_file.write(scenario_name[:-1] + "  " + str(self.__results["Z"])+"\n")
        W_file = open(output_lead+'WW.csv',filewrite)
        X_file = open(output_lead+'X.csv', filewrite)
        L_file = open(output_lead+'L.csv', filewrite)
        #L_bkt_file = open(output_lead+'L_bkt.csv', filewrite)
        
        #Generate WW file
        for gen_name in gen_names: W_file.write(scenario_name+gen_name+"  1 \n")
        for pv_name_panel in pv_name_panels: 
            X_file.write(scenario_name+pv_name_panel+"\n")
        for bat_name in bat_names: W_file.write(scenario_name+bat_name+"  1 \n")
            
        #Generate L file
        for i, gen_name in enumerate(gen_names): 
            L_file.write(scenario_name+gen_name + "  " + str(L_gen[i]) + "\n")
        for i, bat_name in enumerate(bat_names):
            L_file.write(scenario_name+bat_name + "  " + str(L[i]) + "\n")
        
        for t in range(len(self.__results["P_out_g"])):
            #P_out
            for i,p in enumerate(P_out_g[t]): 
                if p > 0.0: P_out_file.write(scenario_name+str(gen_names[i])+'.TT'+str(t+1)+'  '+str(p)+'\n')
            for i,p in enumerate(P_out_s[t]): 
                if p > 0.0: P_solar_file.write(scenario_name+pv_names[i]+'.TT'+str(t+1)+'  '+str(p)+'\n')
            for i,p in enumerate(P_out_b[t]): 
                if p > 0.0: P_out_file.write(scenario_name+bat_names[i]+'.TT'+str(t+1)+'  '+str(p)+'\n')
            #P_in
            for i,p in enumerate(P_in_b[t]): 
                if p > 0.0: P_in_file.write(scenario_name+bat_names[i]+'.TT'+str(t+1)+'  '+str(p)+'\n')
            #B_soc
            for i,p in enumerate(B_soc[t]): 
                if p > 0.0: B_soc_file.write(scenario_name+bat_names[i]+'.TT'+str(t+1)+'  '+str(p)+'\n')
            #B_plus
            for i,p in enumerate(B_plus[t]): 
                if p > 0.0: B_plus_file.write(scenario_name+bat_names[i]+'.TT'+str(t+1)+'  '+str(p)+'\n')
            #B_minus
            for i,p in enumerate(B_minus[t]): 
                if p > 0.0: B_minus_file.write(scenario_name+bat_names[i]+'.TT'+str(t+1)+'  '+str(p)+'\n')
            """
            #L_bkt
            for i,p in enumerate(L_bkt[t]): 
                L_bkt_file.write(scenario_name+bat_names[i]+'.TT'+str(t+1)+'  '+str(p)+'\n')            
            #Y_plus
            for i,p in enumerate(Y_plus[t]): 
                Y_plus_file.write(scenario_name+bat_names[i]+'.TT'+str(t+1)+'  '+str(p)+'\n')
            Y_minus
            for i,p in enumerate(Y_minus[t]): 
                Y_minus_file.write(scenario_name+bat_names[i]+'.TT'+str(t+1)+'  '+str(p)+'\n')
            """
            #Z_plus
            for i,p in enumerate(Z_plus[t]): 
                Z_plus_file.write(scenario_name+bat_names[i]+'.TT'+str(t+1)+'  '+str(p)+'\n')
            #Z_minus
            for i,p in enumerate(Z_minus[t]): 
                Z_minus_file.write(scenario_name+bat_names[i]+'.TT'+str(t+1)+'  '+str(p)+'\n')
            #G
            for i,p in enumerate(G[t]): 
                if p > 0: G_file.write(scenario_name+gen_names[i]+'.TT'+str(t+1)+'  '+str(p)+'\n')
            #G_start
            for i,p in enumerate(G_start[t]): 
                if p > 0: G_start_file.write(scenario_name+gen_names[i]+'.TT'+str(t+1)+'  '+str(p)+'\n')
            #G_stop
            for i,p in enumerate(G_stop[t]): 
                if p > 0: G_stop_file.write(scenario_name+gen_names[i]+'.TT'+str(t+1)+'  '+str(p)+'\n')
            #V_soc
            for i,p in enumerate(V_soc[t]): 
                if p > 0.0: V_soc_file.write(scenario_name+bat_names[i]+'.TT'+str(t+1)+'  '+str(p)+'\n')
            #I_plus
            for i,p in enumerate(I_plus[t]): 
                if p > 0.0: I_plus_file.write(scenario_name+bat_names[i]+'.TT'+str(t+1)+'  '+str(p)+'\n')
            #I_minus
            for i,p in enumerate(I_minus[t]): 
                if p > 0.0: I_minus_file.write(scenario_name+bat_names[i]+'.TT'+str(t+1)+'  '+str(p)+'\n')
            """
            #C_adj
            for i,p in enumerate(C_adj[t]): 
                C_adj_file.write(scenario_name+bat_names[i]+'.TT'+str(t+1)+'  '+str(p)+'\n')
            """
            #F_tilde
            if F_tilde[t] > 0.0: F_tilde_file.write(scenario_name+'TT'+str(t+1)+'  '+str(F_tilde[t])+'\n')
        
    def GenerateAMPLFile(self,scenario_name = None,output_lead = None,append=False):
        if append: filewrite = 'a'
        else: filewrite = 'w'
        #generate additional battery and generator parameters
        battery_binaries = self.CalculateBatteryBinaries()
        #V_soc, I_minus, I_plus = self.CalculateVoltCurrentCap(battery_binaries)
        G_start, G_stop = self.CalculateGeneratorBinaries()
        gen_names = [gen.GetName()+' '+str(gen.GetTwin()) for gen in self.__dispatcher.GetDieselGenerators()]
        pv_names = [pv.GetName() for pv in self.__dispatcher.GetPVArrays()]
        pv_name_panels = [pv.GetName()+'   '+str(pv.GetNumberOfPanels()) for pv in self.__dispatcher.GetPVArrays()]
        bat_names = [bat.GetName() + ' '+str(bat.GetTwin()) for bat in self.__dispatcher.GetBatteries()]
        #Set up variables
        B_soc = self.__results["B_soc"]
        P_out_b = self.__results["P_out_b"]
        P_in_b = self.__results["P_in_b"]
        P_out_s = self.__results["P_out_s"]
        P_out_g = self.__results["P_out_g"]
        G = self.__results["G"]
        L_gen = [sum([G[elt][gg] for elt in range(len(G))]) for gg in range(len(G[0]))]
        L = self.__results["L"]
        #L_bkt = self.__results["L_bkt"]
        F_tilde = self.__results["F_tilde"]
        B_plus = battery_binaries["B_plus"]
        B_minus = battery_binaries["B_minus"]
        I_plus = self.__results["I_plus"]
        I_minus = self.__results["I_minus"]
        Y_plus = self.__results["I_plus"]
        Y_minus = self.__results["I_minus"]
        Z_plus = self.__results["Z_plus"]
        Z_minus = self.__results["Z_minus"]        
        #Set up output file
        outfile = open(output_lead+scenario_name+'_IFS.dat',filewrite) 
        #Generate objective value Z0
        outfile.write("param Obj0 := " + str(self.__results["Z"]))
        outfile.write(";\n\n")
        #Generate W0 values
        outfile.write("param W0 := ")
        for gen_name in gen_names: outfile.write("\n"+gen_name+"  1")
        for bat_name in bat_names: outfile.write("\n"+bat_name+"  1")
        outfile.write(";\n\n")
        #Generate X values
        if len(pv_name_panels) > 0: outfile.write("param X0 := ")
        for pv_name_panel in pv_name_panels: 
            outfile.write("\n"+pv_name_panel)   
        outfile.write(";\n\n") 
        #Generate L values
        outfile.write("param L0 := ")
        for i, gen_name in enumerate(gen_names): 
            outfile.write("\n"+gen_name + "  " + str(L_gen[i]))
        for i, bat_name in enumerate(bat_names):
            outfile.write("\n"+bat_name + "  " + str(L[i]))
        outfile.write(";\n\n") 
        
        #P_out
        times = range(len(self.__results["P_out_g"]))
        outfile.write("param P_minus0 := ")
        for t in times:    
            for i,p in enumerate(P_out_g[t]): 
                if p > 0.0: outfile.write('\n'+str(gen_names[i])+' '+str(t+1)+'  '+str(p))
            for i,p in enumerate(P_out_b[t]): 
                if p > 0.0: outfile.write('\n'+bat_names[i]+' '+str(t+1)+'  '+str(p))
        outfile.write(";\n\n") 
        #P_solar
        outfile.write("param P_PV0 := ")
        for t in times:    
            for i,p in enumerate(P_out_s[t]): 
                if p > 0.0: outfile.write('\n'+pv_names[i]+' '+str(t+1)+'  '+str(p))
        outfile.write(";\n\n") 
        #P_in
        outfile.write("param P_plus0 := ")
        for t in times: 
            for i,p in enumerate(P_in_b[t]): 
                if p > 0.0: outfile.write('\n'+bat_names[i]+' '+str(t+1)+'  '+str(p))
        outfile.write(";\n\n") 
        #B_soc
        outfile.write("param B_soc0 := ")
        for t in times: 
            for i,p in enumerate(B_soc[t]): 
                if p > 0.0: outfile.write('\n'+bat_names[i]+' '+str(t+1)+'  '+str(p))
        outfile.write(";\n\n") 
        #B_plus
        outfile.write("param B_plus0 := ")
        for t in times:
            for i,p in enumerate(B_plus[t]): 
                if p > 0.0: outfile.write('\n'+bat_names[i]+' '+str(t+1)+'  '+str(p))
        outfile.write(";\n\n") 
        #B_minus
        outfile.write("param B_minus0 := ")
        for t in times:
            for i,p in enumerate(B_minus[t]): 
                if p > 0.0: outfile.write('\n'+bat_names[i]+' '+str(t+1)+'  '+str(p))
        outfile.write(";\n\n")
        #Y_minus        
        outfile.write("param Y_minus0 := ")
        for t in times:
            for i,p in enumerate(Y_minus[t]): 
                if p > 0.0: outfile.write('\n'+bat_names[i]+' '+str(t+1)+'  '+str(p))
        outfile.write(";\n\n")
        #Y_minus
        outfile.write("param Y_plus0 := ")
        for t in times:
            for i,p in enumerate(Y_plus[t]): 
                if p > 0.0: outfile.write('\n'+bat_names[i]+' '+str(t+1)+'  '+str(p))
        outfile.write(";\n\n")
        #Z_minus
        outfile.write("param Z_minus0 := ")
        for t in times:
            for i,p in enumerate(Z_minus[t]): 
                if p > 0.0: outfile.write('\n'+bat_names[i]+' '+str(t+1)+'  '+str(p))
        outfile.write(";\n\n")
        #Z_plus
        outfile.write("param Z_plus0 := ")
        for t in times:
            for i,p in enumerate(Z_plus[t]): 
                if p > 0.0: outfile.write('\n'+bat_names[i]+' '+str(t+1)+'  '+str(p))
        outfile.write(";\n\n")
        #G
        outfile.write("param Gen0 := ")
        for t in times:
            for i,p in enumerate(G[t]): 
                if p > 0: outfile.write('\n'+str(gen_names[i])+' '+str(t+1)+'  '+str(p))
        outfile.write(";\n\n")
        """
        #G_start
        outfile.write("param G_start0 := ")
        for t in times:
            for i,p in enumerate(G_start[t]): 
                if p > 0: outfile.write('\n'+str(gen_names[i])+' '+str(t+1)+'  '+str(p))
        #G_stop
        outfile.write("param G_stop0 := ")
        for t in times:
            for i,p in enumerate(G_stop[t]): 
                if p > 0: outfile.write('\n'+str(gen_names[i])+' '+str(t+1)+'  '+str(p))
        #V_soc
            for i,p in enumerate(V_soc[t]): 
                if p > 0.0: V_soc_file.write(scenario_name+bat_names[i]+' '+str(t+1)+'  '+str(p)+'\n')
        """
        #I_plus
        outfile.write("param I_plus0 := ")
        for t in times:
            for i,p in enumerate(I_plus[t]): 
                if p > 0.0: outfile.write('\n'+str(bat_names[i])+' '+str(t+1)+'  '+str(p))
        outfile.write(";\n\n")
        #I_minus
        outfile.write("param I_minus0 := ")
        for t in times:
            for i,p in enumerate(I_minus[t]): 
                if p > 0: outfile.write('\n'+str(bat_names[i])+' '+str(t+1)+'  '+str(p))
        outfile.write(";\n\n")
        #F_tilde
        outfile.write("param F_tilde0 := ")
        for t in times:
            if F_tilde[t] > 0.0: outfile.write('\n'+str(t+1)+'  '+str(F_tilde[t]))
        outfile.write(";\n\n")   
        
        
        
    
        
            

    
