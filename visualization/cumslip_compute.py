#!/usr/bin/env python3
'''
Functions related to plotting cumulative slip vs. depth plot
By Jeena Yun
Last modification: 2023.07.05.
'''
import numpy as np
from scipy import interpolate

yr2sec = 365*24*60*60
wk2sec = 7*24*60*60

def event_times(dep,outputs,Vlb=0,Vths=1e-2,cuttime=0,dt_coseismic=0.5,print_on=True):
    time = np.array(outputs[0][:,0])
    sliprate = abs(np.array([outputs[i][:,4] for i in np.argsort(abs(dep))]))
    z = np.sort(abs(dep))

    if abs(cuttime) >= 1e-3:
        if cuttime > np.max(time):
            raise ValueError('Cuttime larger than total simulation time - check again')
        sliprate = sliprate[:,time <= cuttime]
        time = time[time <= cuttime]

    psr = np.max(sliprate,axis=0)
    pd = np.argmax(sliprate,axis=0)

    # Define events by peak sliprate
    if Vlb > 0:
        events = np.where(np.logical_and(psr < Vths,psr > Vlb))[0]
    else:
        events = np.where(psr > Vths)[0]

    jumps = np.where(np.diff(events)>1)[0]+1

    tmp_tstart = time[events][np.hstack(([0],jumps))]
    tmp_tend = time[events][np.hstack((jumps-1,len(events)-1))]
    tmp_evdep = pd[events][np.hstack(([0],jumps))]

    ii = np.where(tmp_tend-tmp_tstart>=dt_coseismic)[0]
    tstart = tmp_tstart[ii]
    tend = tmp_tend[ii]
    amax = tmp_evdep[ii]
    evdep = z[tmp_evdep[ii]]

    its_all = np.array([np.argmin(abs(time-t)) for t in tstart])
    ite_all = np.array([np.argmin(abs(time-t)) for t in tend])
    diffcrit = np.quantile(abs(np.diff(np.log10(psr[events]))),0.98)
    new_its_all = its_all.copy()
    for k,ts in enumerate(its_all):
        psr_inc = abs(np.diff(np.log10(psr)))[ts-1]
        width = int((ite_all[k] - ts)*0.2)
        large_diffs = np.where(abs(np.diff(np.log10(psr)))[ts-1:ts-1+width]>=diffcrit)[0]
        if psr_inc < diffcrit and len(large_diffs) > 0:
            new_its_all[k] = large_diffs[0] + ts
    evdep = z[pd[new_its_all]]
    tstart = time[new_its_all]

    return tstart, tend, evdep

def compute_cumslip(outputs,dep,cuttime,Vlb,Vths,dt_creep,dt_coseismic,dt_interm,print_on=True):
    if print_on: print('Cumulative slip vs. Depth plot >>> ',end='')

    if print_on: 
        if abs(cuttime) < 1e-3:
            print('No cutting')
        else:
            print('Cut at %2.1f yr'%(cuttime/yr2sec))
        if Vlb > 0:
            print('%1.0e < Slip rate < %1.0e'%(Vlb,Vths))
        else:
            print('Slip rate > %1.0e'%(Vths))

    cscreep = []
    depcreep = []
    cscoseis = []
    depcoseis = []
    fault_slip = []
    if dt_interm > 0:        
        csinterm = []
        depinterm = []

    # Obtain globally min. event start times and max. event tend times
    if dt_interm > 0:
        tstart_interm, tend_interm, evdep = event_times(dep,outputs,Vlb,Vths,cuttime,dt_coseismic,print_on)
    tstart_coseis, tend_coseis, evdep = event_times(dep,outputs,0,Vths,cuttime,dt_coseismic,print_on)
    evslip = np.zeros(tstart_coseis.shape)

    # Now interpolate the cumulative slip using given event time ranges
    for i in np.argsort(abs(dep)):
        z = abs(dep[i])
        time = np.array(outputs[i])[:,0]
        sliprate = np.array(outputs[i])[:,4]
        cumslip = np.array(outputs[i])[:,2]

        if abs(cuttime) >= 1e-3:
            sliprate = sliprate[time <= cuttime]
            cumslip = cumslip[time <= cuttime]
            time = time[time <= cuttime]

        f = interpolate.interp1d(time,cumslip)

        # -------------------- Creep
        tcreep = np.arange(time[0],time[-1],dt_creep)
        cscreep.append(f(tcreep))
        depcreep.append(z*np.ones(len(tcreep)))
        
        # -------------------- Inter
        if dt_interm > 0:
            cs = []
            depth = []
            for j in range(len(tstart_interm)):
                tinterm = np.arange(tstart_interm[j],tend_interm[j],dt_interm)
                cs.append(f(tinterm))
                depth.append(z*np.ones(len(tinterm)))

            csinterm.append([item for sublist in cs for item in sublist])
            depinterm.append([item for sublist in depth for item in sublist])

        # -------------------- Coseismic
        cs = []
        depth = []
        Dbar = []
        for j in range(len(tstart_coseis)):
            tcoseis = np.arange(tstart_coseis[j],tend_coseis[j],dt_coseismic)
            cs.append(f(tcoseis))
            depth.append(z*np.ones(len(tcoseis)))
            Dbar.append(f(tcoseis)[-1]-f(tcoseis)[0])

        cscoseis.append([item for sublist in cs for item in sublist])
        depcoseis.append([item for sublist in depth for item in sublist])
        fault_slip.append(Dbar)

        # -------------------- Event detph
        if np.isin(z,evdep):
            indx = np.where(z==evdep)[0]
            evslip[indx] = f(tstart_coseis)[indx]

    timeout = [tstart_coseis,tend_coseis]
    evout = [evslip,evdep,fault_slip]
    creepout = [cscreep,depcreep]
    coseisout = [cscoseis,depcoseis]
    if dt_interm > 0:
        intermout = [csinterm,depinterm]
    
    if dt_interm > 0:
        return [timeout, evout, creepout, coseisout, intermout]
    else:
        return [timeout, evout, creepout, coseisout]
    
def compute_spinup(outputs,dep,cuttime,cumslip_outputs,spin_up,print_on=True):
    # cumslip_outputs = [timeout, evout, creepout, coseisout, intermout]
    # cumslip_outputs[0] = [tstart_coseis,tend_coseis]
    # cumslip_outputs[1] = [evslip,evdep,fault_slip]
    # cumslip_outputs[2] = [cscreep,depcreep]
    # cumslip_outputs[3] = [cscoseis,depcoseis]
    # cumslip_outputs[4] = [csinterm,depinterm]
    if len(cumslip_outputs) > 4:
        interm = True
    else:
        interm = False
    if print_on: print('Spin-up applied after slip > %2.2f m'%(spin_up))
    spin_up_idx = np.where(cumslip_outputs[1][0]>spin_up)[0][0]
    spup_cscreep = np.copy(cumslip_outputs[2][0])
    if interm:
        spup_csinterm = np.copy(cumslip_outputs[4][0])
    spup_cscoseis = np.copy(cumslip_outputs[3][0])
    spup_evslip = np.copy(cumslip_outputs[1][0])

    new_init_Sl = []
    new_init_dp = []
    c = 0
    for i in np.argsort(abs(dep)):
        z = abs(dep[i])
        time = np.array(outputs[i])[:,0]
        cumslip = np.array(outputs[i])[:,2]

        if abs(cuttime) >= 1e-3:
            sliprate = sliprate[time <= cuttime]
            time = time[time <= cuttime]

        f = interpolate.interp1d(time,cumslip)
        new_init_Sl.append(f(cumslip_outputs[0][0][spin_up_idx]-1))
        new_init_dp.append(z)
        
        spup_cscreep[c] = cumslip_outputs[2][0][c] - new_init_Sl[-1]
        if interm:
            spup_csinterm[c] = cumslip_outputs[4][0][c] - new_init_Sl[-1]
        spup_cscoseis[c] = cumslip_outputs[3][0][c] - new_init_Sl[-1]

        if np.isin(z,cumslip_outputs[1][1]):
            indx = np.where(z==cumslip_outputs[1][1])[0]
            spup_evslip[indx] = cumslip_outputs[1][0][indx] - new_init_Sl[-1]
        c += 1
    
    new_inits = [new_init_Sl,new_init_dp]

    if interm:
        return [new_inits, spup_evslip, spup_cscreep, spup_cscoseis, spup_csinterm, spin_up_idx]
    else:
        return [new_inits, spup_evslip, spup_cscreep, spup_cscoseis, spin_up_idx]
    
def analyze_events(cumslip_outputs,rths):
    from scipy import integrate
    rupture_length = []
    av_slip = []
    fault_z = np.array(cumslip_outputs[3][1]).T[0]
    fault_slip = np.array(cumslip_outputs[1][2]).T

    for ti in range(fault_slip.shape[0]):
        fs = fault_slip[ti]
        Sths = max(fs)*0.01
        rl = max(fault_z[fs>Sths])-min(fault_z[fs>Sths])        
        Dbar = integrate.simpson(fs[fs>Sths],fault_z[fs>Sths])/rl
        rupture_length.append(rl)
        av_slip.append(Dbar)

    rupture_length = np.array(rupture_length)
    partial_rupture = np.where(rupture_length<rths)[0]
    system_wide = np.where(rupture_length>=rths)[0]
    return rupture_length,av_slip,system_wide,partial_rupture