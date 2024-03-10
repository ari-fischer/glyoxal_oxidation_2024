def collect_E(df):
#takes in a df, and then finds the important energy information contained in the paths
    import numpy as np
    #electronic energy
    Es = np.array([])
    #free energy for cavitation
    G_cds = np.array([]) 
    for i in range(df.shape[0]):
        #local path to geometries
        lp = df['path_l_geom'][i]
        print(lp)
        #import energy outputs
        with open(lp+'/Energy.txt') as f:
            lines = f.readlines() 
        if np.size(lines) > 0:
            E = float(lines[0].split()[3])
        else:
            E = 0
        #import cavitation, dispersion, solvation  # 1 Molar reference state 
        with open(lp+'/G_cds.txt') as f:
            lines = f.readlines() 
        if np.size(lines) > 0:
            #print(lines[0].split())
            G = float(lines[0].split()[2])
        else:
            G = 0
        Es = np.append(Es,E)
        G_cds = np.append(G_cds,G)
    
    #add to the dataframe
    df['G_cds'] = G_cds*4.1844 # 1 molar reference state

    df['E'] = Es*2625.5
    return df

    def collect_freq(df,T_C):

#takes in dataframe from CSV with all the important info
#updated to export # of imaginary and the value of the first one

#import information from frequency calculations
# This block calculates the thermodynamic information based on vibrational freuencies as harmonic oscillators
    import numpy as np
    import pandas as pd
    
    #important fundamental constants
    [N_A,k_B,T,h,kJmol_Ha] = [6.0221408e23,1.380649e-23, 273+T_C,6.626070E-34,2625.5] #[ -- , m2 kg s-2 K-1, K,  # m2 kg / s]
    R_ig = 8.3144E-3
    R = k_B*N_A # ideal gas constant: J / mol / K

    h_ev = 4.135668E-15 #planks constant in ev*s
    k_ev = 0.0000861733 #boltzmann in ev/K

    #initial vector to build S H and ZPEs
    S_vec = np.array([])
    H_vec = np.array([])
    ZPE_vec = np.array([])
    H_rot_vec = np.array([])
    H_trans_vec = np.array([])
    S_rot_ig_vec = np.array([])
    S_trans_vec = np.array([])
    S_elec_vec = np.array([])
    n_ims = np.array([])
    vibs_list = np.array([])
    m_red = np.array([])

    #loop over each structure
    for i in range(df.shape[0]):
        
        lp = df['path_l_freq'][i]

        print(lp)
        #open vibrational frequencies
        with open(lp+'/Frequencies.txt') as f:
            lines = f.readlines() 
        #checks that the freq completed and gave answer output
        if np.size(lines) > 0: 
            vibs = pd.read_csv(lp+'/Frequencies.txt',skiprows=0,header=None,delim_whitespace=True)
            #lowest mode
            vib_0 = vibs[0][0]
            #checks if there are imaginary modes
            if df['Calc'][i][0:2]=='TS':
                #for TS excluding one of the modes corresponding to the direction along the rxn coordiante. TS needs to be in the folder name
                n_im = np.sum(vibs[1]<0)
            else:
                n_im = np.sum(vibs[0]<0)
            n_ims = np.append(n_ims,n_im)
                
            #moments of inertia
            I_iner = pd.read_csv(lp+'/inertia.txt',skiprows=0,header=None,delim_whitespace=True).to_numpy()
            
            #apply threshold for low wavenumber vib modes
            vibs[vibs<100]=100 #in cm01

            #convert to hz
            v_hz = vibs.to_numpy()*29979.2458*1e6#cm-1 -> Mhz -> hz

            #calculate entropy of vibrational modes
            S_vibs =R*((h_ev*v_hz/(k_ev*T*(np.exp(h_ev*v_hz/k_ev/T)-1))-np.log(1-np.exp(-h_ev*v_hz/k_ev/T))))

            #calculate ZPE and vib enthalpy
            ZPE = 1/2*h_ev*v_hz
            H_vib = h_ev*v_hz*np.exp(-h_ev*v_hz/k_ev/T)/(1-np.exp(-h_ev*v_hz/k_ev/T))

            #populate table with S H and ZPE values
            if df['Calc'][i][0:2]=='TS': #TS
                #for TS excluding one of the modes corresponding to the direction along the rxn coordiante. TS needs to be in the folder name
                S_vec = np.append(S_vec,np.sum(np.sum(S_vibs))-S_vibs[0][0]) 
                H_vec = np.append(H_vec,np.sum(np.sum(H_vib))-H_vib[0][0]) 
                ZPE_vec = np.append(ZPE_vec,np.sum(np.sum(ZPE))-ZPE[0][0]) 
                
            else: #minima
                #for minima keep all modes
                S_vec = np.append(S_vec,np.sum(np.sum(S_vibs))) 
                H_vec = np.append(H_vec,np.sum(np.sum(H_vib))) 
                ZPE_vec = np.append(ZPE_vec,np.sum(np.sum(ZPE))) 
            
            # H rot and trans from QChem
            df_H = pd.read_csv(lp + '/H_out'  +'.txt',skiprows=0,header=None,delim_whitespace=True)*4.184
            # S trans, rot, vib, tot from QChem
            df_S = pd.read_csv(lp + '/S_out'  +'.txt',skiprows=0,header=None,delim_whitespace=True)*4.184
            
            #populate the table
            #https://cccbdb.nist.gov/thermox.asp
            #for non-lilnear molecules (need to change for OH and O2)
            H_rot = 3/2*R*T/1000 # kJ/mol 
            #for polyatomic gas (diff for H)
            H_trans = 5/2*R*T/1000#df_H.to_numpy()[0] # kJ/mol

            #from QChem output
            S_rot_ig = df_S.to_numpy()[1] # J/mol/K     
            S_trans = df_S.to_numpy()[0] # J/mol/K   

            #Get electronic entropy from spin multiplicity. Read out_head (with charge and spin info)
            with open(df['path_l_geom'][i]+'/out_head.txt') as f:
                lines = f.readlines()
            M = int(lines[-1][-2])        
        
            if df['Calc'][i]=='OH':
                S_elec = R*np.log(4)
            else:
                S_elec = R*np.log(M)

            #append vectors
            H_rot_vec = np.append(H_rot_vec,H_rot)
            H_trans_vec = np.append(H_trans_vec,H_trans)
            S_rot_ig_vec = np.append(S_rot_ig_vec,S_rot_ig)
            S_trans_vec = np.append(S_trans_vec,S_trans)
            S_elec_vec = np.append(S_elec_vec,S_elec)
            
        else: #populates with zero if there is no computed frequencies
            S_vec = np.append(S_vec,0) 
            #Data['H_vib'][i] = np.sum(np.sum(S_damp))
            H_vec = np.append(H_vec,0) 
            ZPE_vec = np.append(ZPE_vec,0)
            
            H_rot_vec = np.append(H_rot_vec,0)
            H_trans_vec = np.append(H_trans_vec,0)
            S_rot_ig_vec = np.append(S_rot_ig_vec,0)
            S_trans_vec = np.append(S_trans_vec,0)
            S_elec_vec = np.append(S_elec_vec,0)
            n_ims = np.append(n_ims,[0]) 

        vibs_list = np.append(vibs_list,vib_0)
            

    #populate the dataframe
    df['vib_0'] = vibs_list   
    df['n_imag'] = n_ims
    df['S_vib'] = S_vec
    df['H_vib'] = H_vec*96.4869
    df['ZPE']= ZPE_vec*96.4869 #kJ/mol

    df['H_rot'] = H_rot_vec #kcal/mol -> kJ/mol (occured at import)
    df['H_trans'] = H_trans_vec #kcal/mol -> kJ/mol
    df['S_rot'] = S_rot_ig_vec #kcal/mol -> J/mol/K
    # convert from 1 atm to 23.45309 atm (for 1 molar reference state)
    df['S_trans'] = S_trans_vec-np.log(23.4530933)*R_ig*1000 #kcal/mol -> J/mol/K 
    df['S_elec'] = S_elec_vec #kcal/mol -> J/mol/K

    return df

def species_list(s):
    s = [item for sublist in s for item in sublist]

    species = []
    for item in s:
        if '.' in item:
            species.extend(item.split('.'))
        else:
            species.append(item)
        
    species = list(set(species))
    #for item in species:
    #    print(item)
    return species

def rxn_strings(rcts,prod):
    s_rxns = list([])
    for i in range(len(rcts)):
        if len(rcts[i])==2:
            s_rcts = rcts[i][0]+'+'+rcts[i][1]
        elif len(rcts[i])==3:
            s_rcts = rcts[i][0]+'+'+rcts[i][1]+'+'+rcts[i][2]
        else:
            s_rcts = rcts[i][0]
        
        if len(prod[i])==2:
            s_prd = prod[i][0]+'+'+prod[i][1]
        elif len(prod[i])==3:
            s_prd = prod[i][0]+'+'+prod[i][1]+'+'+prod[i][2]
        elif len(prod[i])==4:
            s_prd = prod[i][0]+'+'+prod[i][1]+'+'+prod[i][2]+'+'+prod[i][3]    
        else:
            s_prd = prod[i][0]
        
        #print(s_rcts+'->'+s_prd)
        s_rxns.append(s_rcts+'->'+s_prd)
    return s_rxns
