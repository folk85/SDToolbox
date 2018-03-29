'''
 Shock and Detonation Toolbox
 http://www.galcit.caltech.edu/EDL/public/cantera/html/SD_Toolbox/
'''

from SDToolbox import *
from .CJ2 import *
import datetime

class cvoutput:
    def __init__(self,size=1,nsp=2):
        self.exo_time = 0
        self.ind_time = 0
        self.ind_time_10 = 0
        self.ind_time_90 = 0
        
        self.time = zeros(size,float)
        self.T = zeros(size,float)
        self.P = zeros(size,float)
        self.species = zeros((size,nsp),float)

def cv_CJ(plt_num, P1, T1, q, mech, fname):
    """

    cv_CJ.m
    Computes the time evolution of a constant volume explosion with
    shock (at CJ speed) heated reactants as initial conditions
    
    FUNCTION
    SYNTAX
    [gas] = cv_CJ(plt_num,P1,T1,q,mech,fname)
    
    INPUT
    plt_num = unused
    P1 = initial pressure (Pa)
    T1 = initial temperature (K)
    q = string of reactant species mole fractions
    mech = cti file containing mechanism data (i.e. 'gri30.cti')
    fname = output file name (0 for no output file)
    
    OUTPUT
    gas = gas object at final equilibrium state

    """
    gas1 = Solution(mech);
    gas1.TPX = T1,P1,q; 
    r = gas1.density;
    [gas, cj_speed] = CJspeed2(P1, T1, q, mech);
    gas = PostShock_fr(cj_speed, P1, T1, q, mech);
    fname = fname + '_%d' % cj_speed
    b = 10000; j = gas.n_species;
    out = cvoutput(b,j)
    out = explosion(gas,fname,out);
    return [cj_speed, gas]

def cv_shk(U1, P1, T1, q, mech, fname):
    """

    cv_shk.m
    Computes the time evolution of a constant volume explosion with
    shock heated reactants as initial conditions
    
    FUNCTION
    SYNTAX
    [gas] = cv_shk(U1,P1,T1,q,mech,fname)
    
    INPUT
    U1 = shock speed (m/s)
    P1 = initial pressure (Pa)
    T1 = initial temperature (K)
    q = string of reactant species mole fractions
    mech = cti file containing mechanism data (i.e. 'gri30.cti')
    fname = output file name (0 for no output file)
    
    OUTPUT
    gas = gas object at final equilibrium state

    """
    gas1 = Solution(mech);
    gas1.TPX = T1,P1,q; 
    r = gas1.density;
    gas = PostShock_fr(U1, P1, T1, q, mech);
    fname = fname + '_%d' % U1
    b = 10000; j = gas1.n_species;
    out = cvoutput(b,j)
    out = explosion(gas,fname,out);
    return gas

def explosion(gas,fname,out):
    """

    explosion.m
    Computes the time evolution of a constant volume explosion
    
    FUNCTION
    SYNTAX
    [ind_time,ind_time_10,ind_time_90,time,temp,press,species] = explosion(gas,fig_num)
    
    INPUT
    gas = working gas object
    fig_num = figure number (0 for no plot)
    
    OUTPUT
    exo_time = pulse width (in secs) of temperature gradient (using 1/2 max)
    ind_time = time to maximum temperature gradient
    ind_len = distance to maximum temperature gradient
    ind_time_10 = time to 10% of maximum temperature gradient
    ind_time_90 = time to 90% of maximum temperature gradient
    time = trray of time
    temp = temperature profile array
    press = pressure profile array
    species = matrix of species profiles

    """
    b = 10000; j = gas.n_species; rho = gas.density
    r = Reactor(gas)
    sim = ReactorNet([r])
    t = 0.0
    temp_grad = zeros(b,float)
    y = zeros(j,float)

    #EXTRACT OUTPUT INFORMATION
    for n in range(b):
        out.time[n] = t
        out.T[n] = r.thermo.T
        out.P[n] = r.thermo.P
        for i in range(j):
            out.species[n,i] = r.Y[i]
            y[i] = r.Y[i]
        
        gas.TDY = out.T[n],r.density,y
        P = gas.P/one_atm

        # FIND TEMPERATURE GRADIENT
        # Conservation of Energy implies that e0 = e(T)
        # e = cv*T; dedt = 0; cv*dTdt + sum(ei(T)*dyidt) = 0; dyidt = wdoti*wti/rho
        # dTdt = -sum(ei(T)*wdoti*wti)/(cv*rho)
        cv = gas.cv_mass;
        wdot = gas.net_production_rates;
        mw = gas.molecular_weights
        hs = gas.standard_enthalpies_RT
        R = gas_constant;
        wt = gas.mean_molecular_weight
        
        sumT = 0.0
        for z in range(j):
            w = mw[z]; e = R*out.T[n]*(hs[z]/w - 1/wt)
            wd = wdot[z]; sumT = sumT + e*wd*w;
        temp_grad[n] = -sumT/(rho*cv)
        t += 1.e-9
        sim.advance(t)

    del sim
    del r
    
    #FIND INDUCTION TIME - MAXIMUM TEMPERATURE GRADIENT
    k = 0; MAX = max(temp_grad); d = temp_grad[0]; HMPWt = zeros(2,float)
    if d == MAX:
        print('Initial Gradient is Maximum - post shock temperature may be too low')
        return gas
    while d < MAX:
        k = k + 1; d = temp_grad[k];
    out.ind_time = out.time[k]; k1 = k; k = 0;
    MAX10 = 0.1*MAX; d = temp_grad[0];
    while(d < MAX10 and k < b-1):
        k = k + 1; d = temp_grad[k];
    if(k == b):
        print('MAX10 may be incorrect - reached end of array')
    out.ind_time_10 = out.time[k]; k = 0;
    MAX90 = 0.9*MAX; d = temp_grad[0];
    while(d < MAX90 and k < b-1):
        k = k + 1; d = temp_grad[k];
    if(k == b-1):
        print('MAX90 may be incorrect - reached end of array')
    out.ind_time_90 = out.time[k];

    # find exothermic time
    half_T_flag1 = 0;
    half_T_flag2 = 0;
    #Go into a loop to find two times when Temperature is half its maximum
    for j in range(b) : 
        if (half_T_flag1 == 0 ):
            if (temp_grad[j] >= (0.5* MAX)):
                half_T_flag1 = 1;
                tstep1 = j;
        else:
            if (half_T_flag2 == 0) :
                if (temp_grad[j] <= (0.5* MAX) ):
                    half_T_flag2 = 1;
                    tstep2 = j;
                else:
                    tstep2 = 0;
    #Exothermic time for constant volume explosion
    out.exo_time = out.time[tstep2] - out.time[tstep1];  
    
    if fname==0:
        return out

    else:
        k = 0; MAX = max(out.T); d = out.T[0];
        while d < MAX:
            k = k + 1; d = out.T[k];
        if out.time[k] == 0:
            maxt = out.ind_time*5;
        elif out.time[k] >= out.ind_time*50:
            maxt = out.ind_time*5;
        else:
            maxt = out.time[k] + 0.1*out.time[k];
            mint = 0;
        maxT = max(out.T)+0.1*min(out.T); minT = min(out.T)-0.1*min(out.T); 
        maxP = max(out.P)+0.1*min(out.P); minP = min(out.P)-0.1*min(out.P); 
        maxpw = HMPWt[1] + 0.1*HMPWt[1]; minpw = HMPWt[0] - 0.1*HMPWt[0]; 
        maxTG = max(temp_grad) + 0.1*abs(max(temp_grad));
        minTG = min(temp_grad)-0.1*abs(min(temp_grad));
        d = datetime.date.today(); P = out.P[0]/one_atm;

        fn = fname + '_CVprofile.plt';
        outfile = file(fn, 'w');
        outfile.write('# CONSTANT VOLUME PROFILES\n');
        outfile.write('# CALCULATION RUN ON %s\n\n' % d);
        outfile.write('# Maximum time calculated = %.4e\n' % max(out.time))
        outfile.write('# t_min = %.4f, t_max = %.4e\n' % (mint, maxt))
        outfile.write('# T_min = %.2f, T_max = %.2f\n' % (minT, maxT))
        outfile.write('# P_min = %.2f, P_max = %.2f\n' % (minP, maxP))
        outfile.write('# TG_min = %.2f, TG_max = %.2f\n' % (minTG, maxTG))
        outfile.write('# THE OUTPUT DATA COLUMNS ARE:\n');
        outfile.write('Variables = "Time", "Temperature", "Pressure", "temp_grad"\n');
        for i in range(b):
            outfile.write('%.4E \t %.4E \t %.4E \t %.4E\n'% (out.time[i],out.T[i],out.P[i],temp_grad[i]));

    return out
     
