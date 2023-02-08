#!/usr/bin/env python3
'''
An executable plotting script for Tandem to save figures directly from a remote server
By Jeena Yun
Last modification: 2023.02.07.
'''
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import glob
import csv

mypink = (230/255,128/255,128/255)
myblue = (118/255,177/255,230/255)
myburgundy = (214/255,0,0)
mynavy = (17/255,34/255,133/255)
yr2sec = 365*24*60*60
wk2sec = 7*24*60*60

# Set input parameters -------------------------------------------------------------------------------------------------------------------
# ***
# save_dir: directory to output files (and to save plots). You may need to modify load/save directories.
# compute: choose either to process the raw output or load processed output (to save time)
# save_figs: select which figure you want to save
#            For the first three plots, non-zero values would indicate its desired depth (surface should be 0.01)
# plot_in_sec: choose whether you want the time axis to be in year or in seconds
# Vths, Vlb: lower- and upper bound for the slip rate. [m/s] Used when plotting the cumslip plot. 
#            Vub will be used as coseismic slip rate threshold. Give Vlb = 0 when intermediate section is unnecessary.
# dt_creep: Contour line interval for the creeping section. Only required when cumslip plotting is on [s]
# dt_coseismic: Contour line interval for the coseismic section. Only required when cumslip plotting is on [s]
# mingap: minimum seperation time between two different events. Used only for plotting purpose.
# dt_interm (optional): Contour line interval for the intermediate section. Only required when cumslip plotting is on [s]. Give zero if unnecessary.
# cuttime (optional): show result up until to the given time to save computation time. Give zero if unnecessary.
# Wf (optional): depth of fault. Used only for plotting purpose (to set ylim). Give zero if unnecessary.
# ***
save_dir = '/export/dump/jyun/DZ_triangular'
compute = 0
save_figs = [0.01, 2, 0, 0, 1] # Slip rate vs. Time | Slip vs. Time | Stress vs. Time | Initial stress | Cumulative slip vs. Depth
plot_in_sec = 0
Vths = 1e-2
Vlb = 1e-8
dt_creep = 2*yr2sec
dt_coseismic = 1
mingap = 60
dt_interm = 1*wk2sec
cuttime = 700*yr2sec
Wf = 24

# Extract data ---------------------------------------------------------------------------------------------------------------------------
if compute:
    print('Compute on - extract outputs...',end=' ')

    fnames = glob.glob('%s/outputs/*dp*.csv'%(save_dir))
    if len(fnames) == 0:
        raise NameError('No such file found - check the input')
    
    outputs = ()
    dep = []

    for fn in fnames:
        with open(fn, 'r') as csvfile:
            csvreader = csv.reader(csvfile)
            stloc = next(csvreader)
            r_x = float(stloc[0].split('[')[-1])
            r_z = float(stloc[1].split(']')[0])
            dep.append(r_z)
                
            next(csvreader)

            dat = []
            for row in csvreader:
                dat.append(np.asarray(row).astype(float))
        
        outputs = outputs + (dat,)
    dep = np.array(dep)
    print('done!')

    print('Save data...',end=' ')
    np.save('%s/outputs'%(save_dir),outputs)
    np.save('%s/outputs_depthinfo'%(save_dir),dep)
    print('done!')
else:
    print('Load saved data: %s/outputs'%(save_dir))
    print('Load saved data: %s/outputs_depthinfo'%(save_dir))
    outputs = np.load('%s/outputs.npy'%(save_dir))
    dep = np.load('%s/outputs_depthinfo.npy'%(save_dir))

# Slip rate vs. Time ---------------------------------------------------------------------------------------------------------------------
if abs(save_figs[0])>0:
    target_depth = save_figs[0] # in km
    indx = np.argmin(abs(abs(dep) - abs(target_depth)))
    print('Slip rate vs. Time plot >> Depth = %1.1f [km]'%abs(dep[indx]))

    plt.rcParams['font.size'] = '15'
    plt.figure(figsize=(8,6))
    if plot_in_sec:        # --- Plot in seconds
        if np.all(np.array(outputs[indx])[:,4]>0):
            plt.plot(np.array(outputs[indx])[:,0], np.log10(np.array(outputs[indx])[:,4]), color='k', lw=2.5)
        else:
            print('Negative slip rate - taking absolute')
            plt.plot(np.array(outputs[indx])[:,0], np.log10(abs(np.array(outputs[indx])[:,4])), color='k', lw=2.5)
        plt.xlabel('Time [s]',fontsize=17)
    else:        # --- Plot in years
        if np.all(np.array(outputs[indx])[:,4]>0):
            plt.plot(np.array(outputs[indx])[:,0]/yr2sec, np.log10(np.array(outputs[indx])[:,4]), color='k', lw=2.5)
        else:
            print('Negative slip rate - taking absolute')
            plt.plot(np.array(outputs[indx])[:,0]/yr2sec, np.log10(abs(np.array(outputs[indx])[:,4])), color='k', lw=2.5)
        plt.xlabel('Time [yr]',fontsize=17)
    plt.ylabel('log$_{10}$(Slip Rate [m/s])',fontsize=17)
    if target_depth < 1e-1:
        plt.title('Depth = surface',fontsize=20,fontweight = 'bold')
    else:
        plt.title('Depth = %1.1f [km]'%target_depth,fontsize=20,fontweight = 'bold')
    plt.tight_layout()
    plt.savefig('%s/sliprate.png'%(save_dir),dpi=300)

# Slip vs. Time --------------------------------------------------------------------------------------------------------------------------
if abs(save_figs[1])>0:
    target_depth = save_figs[1] # in km
    indx = np.argmin(abs(abs(dep) - abs(target_depth)))
    print('Slip vs. Time plot >> Depth = %1.1f [km]'%abs(dep[indx]))
    time = np.array(outputs[indx])[:,0]
    cumslip = np.array(outputs[indx])[:,2]

    plt.rcParams['font.size'] = '15'
    plt.figure(figsize=(8,6))
    if plot_in_sec:        # --- Plot in seconds
        plt.plot(time,cumslip, color='k', lw=2.5)
        plt.xlabel('Time [s]',fontsize=17)
    else:        # --- Plot in years
        plt.plot(time/yr2sec,cumslip, color='k', lw=2.5)
        plt.xlabel('Time [s]',fontsize=17)
    plt.ylabel('Cumulative Slip [m]',fontsize=17)
    if target_depth < 1e-1:
        plt.title('Depth = surface',fontsize=20,fontweight = 'bold')
    else:
        plt.title('Depth = %1.1f [km]'%target_depth,fontsize=20,fontweight = 'bold')
    plt.tight_layout()
    plt.savefig('%s/slip.png'%(save_dir))

# Stress vs. Time ------------------------------------------------------------------------------------------------------------------------
if abs(save_figs[2])>0:
    target_depth = save_figs[2] # in km
    indx = np.argmin(abs(abs(dep) - abs(target_depth)))
    print('Stress vs. Time plot >> Depth = %1.1f [km]'%abs(dep[indx]))

    plt.rcParams['font.size'] = '15'
    fig,ax = plt.subplots(ncols=2,figsize=(14,6))
    if plot_in_sec:        # --- Plot in seconds
        ax[0].plot(np.array(outputs[indx])[:,0], np.array(outputs[indx])[:,3], color='k', lw=2.5)
        ax[0].set_xlabel('Time [s]',fontsize=17)
        ax[1].plot(np.array(outputs[indx])[:,0], np.array(outputs[indx])[:,5], color='k', lw=2.5)
        ax[1].set_xlabel('Time [s]',fontsize=17)
    else:        # --- Plot in years
        ax[0].plot(np.array(outputs[indx])[:,0]/yr2sec, np.array(outputs[indx])[:,3], color='k', lw=2.5)
        ax[0].set_xlabel('Time [yr]',fontsize=17)
        ax[1].plot(np.array(outputs[indx])[:,0]/yr2sec, np.array(outputs[indx])[:,5], color='k', lw=2.5)
        ax[1].set_xlabel('Time [yr]',fontsize=17)
    ax[0].set_ylabel('Shear Stress [MPa]',fontsize=17)
    ax[1].set_ylabel('Normal Stress [MPa]',fontsize=17)
    if target_depth < 1e-1:
        ax[0].set_title('Depth = surface',fontsize=20,fontweight = 'bold')
        ax[1].set_title('Depth = surface',fontsize=20,fontweight = 'bold')
    else:
        ax[0].set_title('Depth = %1.1f [km]'%target_depth,fontsize=20,fontweight = 'bold')
        ax[1].set_title('Depth = %1.1f [km]'%target_depth,fontsize=20,fontweight = 'bold')
    plt.tight_layout()
    plt.savefig('%s/stresses.png'%(save_dir))

# Initial stress vs. Time ------------------------------------------------------------
if save_figs[3]:
    print('Initial stress vs. Time plot')
    z = np.zeros(len(dep))
    tau = np.zeros(len(dep))
    sigma = np.zeros(len(dep))

    c = 0
    for i in np.argsort(abs(dep)):
        z[c] = abs(dep[i])
        tau[c] = np.array(outputs[i])[0,3]
        sigma[c] = np.array(outputs[i])[0,5]
        c += 1

    plt.rcParams['font.size'] = '15'
    plt.figure(figsize=(9,7))
    plt.scatter(sigma,z,lw=2.5,color=mynavy,label='Normal Stress',zorder=3)
    plt.scatter(tau,z,lw=2.5,color=myburgundy,label='Shear Stress',zorder=3)
    plt.xlabel('Stress [MPa]',fontsize=17)
    plt.ylabel('Depth [km]',fontsize=17)
    if abs(Wf) > 1e-3:
        plt.ylim(0,Wf)
    plt.gca().invert_yaxis()
    plt.grid(True)
    plt.legend(fontsize=13,loc='lower left')
    plt.tight_layout()
    plt.savefig('%s/stressvsdepth.png'%(save_dir))

# Cumslip vs. Depth ----------------------------------------------------------------------------------------------------------------------
if save_figs[4]:
    print('Cumulative slip vs. Depth plot >>> ',end='')
    cscreep = []
    depcreep = []
    cscoseis = []
    depcoseis = []
    if dt_interm > 0:        
        csinterm = []
        depinterm = []

    if abs(cuttime) < 1e-3:
        print('No cutting')
    else:
        print('Cut at %2.1f yr'%(cuttime/yr2sec))

    def event_times(dep,outputs,Vlb,Vths,cuttime):
        '''
        Search for the earliest start of the event and the last stop of the event.
        (EXPLANATION - you may read it once and remove it as you wish)
        My strategy to make this plot is to do 1d interpolation of the given cumulative slip data using
        a evenly-spaced time variable for both creeping and coseismic period.
        Creeping period is usually imaged sparsely (e.g., several years), so it does not cause any problem.
        The key to this function is to 
        1) extract event start/end times for each depth 
        2) come up with an array of event start/end times that embraces the start/end times from all depths
        I am doing step 2 because otherwise, the length of the start/end time varies, which causes problem
        when plotting them.
        '''
        c = 0
        for i in np.argsort(abs(dep)):
            z = abs(dep[i])
            time = np.array(outputs[i])[:,0]
            sliprate = np.array(outputs[i])[:,4]
            cumslip = np.array(outputs[i])[:,2]

            if abs(cuttime) < 1e-3:
                if c == 0:
                    print('No cutting')
            else:
                if cuttime > np.max(time):
                    raise ValueError('Cuttime larger than total simulation time - check again')
                else:
                    if c == 0:
                        print('Cut at %2.1f yr'%(cuttime/yr2sec))
                    sliprate = sliprate[time <= cuttime]
                    cumslip = cumslip[time <= cuttime]
                    time = time[time <= cuttime]

            # Define events by sliprate
            if Vlb > 0:
                if c == 0:
                    print('%1.0e < Slip rate < %1.0e'%(Vlb,Vths))
                events = np.where(np.logical_and(sliprate < Vths,sliprate > Vlb))[0]
            else:
                if c == 0:
                    print('Slip rate > %1.0e'%(Vths))
                events = np.where(sliprate > Vths)[0]

            if len(events) == 0:
                print('Depth',z,' - no events')
                continue
            else:
                # Get indexes for the dynamic rupture components
                jumps = np.where(np.diff(events)>1)[0]+1

                # Get event start/end time for current depth
                tmp_tstart = np.zeros(len(jumps)+1)
                tmp_tend = np.zeros(len(jumps)+1)
                tmp_evdep = z*np.ones(len(jumps)+1)
                if len(jumps) > 0:
                    for j in range(len(jumps)+1):
                        if j==0:
                            tmp_tstart[j] = time[events][0]
                            tmp_tend[j] = time[events][jumps[0]-1]
                        elif j == len(jumps):
                            tmp_tstart[j] = time[events][jumps[j-1]]
                            tmp_tend[j] = time[events][len(events)-1]
                        else:
                            tmp_tstart[j] = time[events][jumps[j-1]]
                            tmp_tend[j] = time[events][jumps[j]-1]
                else:
                    continue
                
                if c == 0:
                    # When first depth, initiate tstart
                    tstart = np.copy(tmp_tstart)
                    tend = np.copy(tmp_tend)
                    evdep = np.copy(tmp_evdep)
                else:
                    if len(tmp_tstart) > len(tstart):
                        # More number of events
                        long_tstart = np.copy(tmp_tstart)
                        long_tend = np.copy(tmp_tend)
                        long_evdep = np.copy(tmp_evdep)
                        short_tstart = np.copy(tstart)
                        short_tend = np.copy(tend)
                        short_evdep = np.copy(evdep)
                    elif len(tmp_tstart) <= len(tstart):
                        # Less number of events
                        long_tstart = np.copy(tstart)
                        long_tend = np.copy(tend)
                        long_evdep = np.copy(evdep)
                        short_tstart = np.copy(tmp_tstart)
                        short_tend = np.copy(tmp_tend)
                        short_evdep = np.copy(tmp_evdep)

                    # Iteratively update current event start/end time
                    new_tstart = np.copy(long_tstart)
                    new_tend = np.copy(long_tend)
                    new_evdep = np.copy(long_evdep)

                    want_append = 0
                    for k in range(len(short_tstart)):
                        # Cases when the event need to be inserted
                        if short_tstart[k] > max(long_tend):
                            want_append = 1
                        else:
                            same_event = np.argmin(abs(short_tstart[k] - long_tstart))
                            if long_tstart[same_event] > short_tstart[k]:
                                if same_event > 0:
                                    if short_tstart[k] > long_tend[same_event-1] and short_tend[k] < long_tstart[same_event]:
                                        want_append = 1
                                    elif short_tstart[k] < long_tend[same_event-1] and short_tend[k] > long_tstart[same_event]:
                                        want_append = 1
                                if same_event == 0 and short_tend[k] < long_tstart[0]:
                                    want_append = 1
                            elif long_tstart[same_event] < short_tstart[k]:
                                if same_event+1 < len(long_tstart):
                                    if short_tstart[k] > long_tend[same_event] and short_tend[k] < long_tstart[same_event+1]:
                                        want_append = 1
                                    elif short_tstart[k] < long_tend[same_event] and short_tend[k] > long_tstart[same_event+1]:
                                        want_append = 1
                                if same_event+1 == len(long_tstart) and short_tstart[k] > long_tend[-1]:
                                    want_append = 1            

                            if want_append:
                                new_tstart = np.append(new_tstart,short_tstart[k])
                                new_tend = np.append(new_tend,short_tend[k])
                                new_evdep = np.append(new_evdep,short_evdep[k])
                            else:
                                # Cases when the either start or end time needs to be adjusted
                                if long_tstart[same_event] > short_tstart[k]:
                                    if same_event == 0:
                                        new_tstart[same_event] = short_tstart[k]
                                        new_evdep[same_event] = short_evdep[k]
                                    elif short_tstart[k] > long_tend[same_event-1]:
                                        new_tstart[same_event] = short_tstart[k]
                                        new_evdep[same_event] = short_evdep[k]
                                        if short_tend[k] > long_tend[same_event]:
                                            new_tend[same_event] = short_tend[k]
                                    elif short_tstart[k] < long_tend[same_event-1] and short_tend[k] > long_tend[same_event-1]:
                                        new_tend[same_event-1] = short_tend[k]

                                elif long_tstart[same_event] < short_tstart[k]:
                                    if short_tstart[k] > long_tend[same_event]:
                                        new_tstart[same_event+1] = short_tstart[k]
                                        new_evdep[same_event+1] = short_evdep[k]
                                        if same_event+1 < len(long_tstart) and short_tend[k] > long_tend[same_event+1]:
                                            new_tend[same_event+1] = short_tend[k]
                                    elif short_tstart[k] < long_tend[same_event] and short_tend[k] > long_tend[same_event]:
                                        new_tend[same_event] = short_tend[k]

                                elif long_tstart[same_event] == short_tstart[k]:
                                    if short_tend[k] > long_tend[same_event]:
                                        new_tend[same_event] = short_tend[k]
                        want_append = 0
                    
                    # Sort and update start/end time
                    ii = np.argsort(new_tstart)
                    tstart = np.copy(new_tstart[ii])
                    tend = np.copy(new_tend[ii])
                    evdep = np.copy(new_evdep[ii])
            c += 1

        # If there are too close events, merge them as one
        new_tstart = []
        new_tend = []
        new_evdep = []
        u = 0
        while u < len(tstart):
            nearest = np.where(tstart[tstart > tstart[u]] - tstart[u] <= mingap)[0] + u + 1
            if len(nearest) != 0:
                new_tstart.append(tstart[u])
                new_tend.append(tend[max(nearest)])
                new_evdep.append(evdep[u])
                u += len(nearest)+1
            else:
                new_tstart.append(tstart[u])
                new_tend.append(tend[u])
                new_evdep.append(evdep[u])
                u += 1

        tstart = np.copy(new_tstart)
        tend = np.copy(new_tend)
        evdep = np.copy(new_evdep)

        return tstart, tend, evdep

    # Obtain globally min. event start times and max. event tend times
    if dt_interm > 0:
        tstart_interm, tend_interm, evdep = event_times(dep,outputs,Vlb,Vths,cuttime)
    tstart_coseis, tend_coseis, evdep = event_times(dep,outputs,0,Vths,cuttime)
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
        for j in range(len(tstart_coseis)):
            tcoseis = np.arange(tstart_coseis[j],tend_coseis[j],dt_coseismic)
            cs.append(f(tcoseis))
            depth.append(z*np.ones(len(tcoseis)))

        cscoseis.append([item for sublist in cs for item in sublist])
        depcoseis.append([item for sublist in depth for item in sublist])

        # -------------------- Event detph
        if np.isin(z,evdep):
            indx = np.where(z==evdep)[0]
            evslip[indx] = f(tstart_coseis)[indx]

    # --- Plot the result
    plt.rcParams['font.size'] = '27'
    plt.figure(figsize=(18,11))

    plt.plot(cscreep,depcreep,color='royalblue',lw=1)
    if dt_interm > 0:
        plt.plot(csinterm,depinterm,color='yellowgreen',lw=1)
    plt.plot(cscoseis,depcoseis,color='chocolate',lw=1)
    ev = plt.scatter(evslip,evdep,marker='*',s=700,facecolor=myburgundy,edgecolors='k',lw=2,zorder=3)
    plt.legend([ev],['Hypocenters'],fontsize=25,framealpha=1)
    plt.ylabel('Depth [km]',fontsize=30)
    plt.xlabel('Cumulative Slip [m]',fontsize=30)
    xl = plt.gca().get_xlim()
    plt.xlim(0,xl[1])
    if abs(Wf) > 1e-3:
        plt.ylim(0,Wf)
    plt.gca().invert_yaxis()

    plt.tight_layout()
    plt.savefig('%s/cumslip_%d_%d.png'%(save_dir,int(Vths*100),int(dt_coseismic*10)),dpi=300)

