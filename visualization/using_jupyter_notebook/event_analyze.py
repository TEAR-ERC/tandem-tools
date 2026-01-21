#!/usr/bin/env python3
'''
Functions related to event analyzation
By Jeena Yun
Last modification: 2025.11.18.
'''
import numpy as np

yr2sec = 365*24*60*60
wk2sec = 7*24*60*60

def detect_events(outputs, dep, options):
    '''
    Detecting seismic events
    tstart, tend, evdep:
    evslip: cumulative slip at the start of the event
    rupture_length
    av_slip
    system_wide
    partial_rupture
    event_cluster
    lead_fs, major_pr, minor_pr
    '''
    options['SSE'] = False
    tstart, tend, evdep, its_all, ite_all, Vmax, duration = event_times(outputs, dep, options)

    fault_slip, evslip, rupture_length, av_slip, system_wide, partial_rupture, event_cluster, lead_fs, major_pr, minor_pr = \
         [],[],[],[],[],[],[],[],[],[]
    if len(tstart) > 0:
        for i in range(len(dep)):
            z = abs(dep)[i]
            cumslip = abs(np.array(outputs[i])[:,2])

            # -------------------- Coseismic
            Dbar = []
            for j in range(len(its_all)):
                Dbar.append(cumslip[ite_all[j]] - cumslip[its_all[j]])
                # -------------------- Event detph
                if z == evdep[j]:
                    evslip.append(cumslip[its_all[j]])
            fault_slip.append(Dbar)
        fault_slip = np.array(fault_slip)
        evslip = np.array(evslip)

        rupture_length, av_slip, system_wide, partial_rupture, event_cluster, lead_fs, major_pr, minor_pr = \
            analyze_events(dep, fault_slip, tstart, options)
        
    event_info = {'tstart':tstart, 'tend':tend, 'its_all':its_all, 'ite_all':ite_all,\
                'evdep':evdep, 'evslip':evslip, 'fault_slip':fault_slip,\
                'duration':duration, 'Vmax':Vmax, \
                'rupture_length':rupture_length, 'av_slip':av_slip, 'system_wide':system_wide, 'partial_rupture':partial_rupture, \
                'event_cluster':event_cluster, 'lead_fs':lead_fs, 'major_pr':major_pr, 'minor_pr':minor_pr, \
                'options':options}
    
    return event_info


def detect_SSEs(outputs, dep, options):
    '''
    Detecting SSEs
    tstart, tend, evdep:
    evslip: cumulative slip at the start of the event
    rupture_length
    av_slip
    duration
    '''
    options['SSE'] = True
    tstart, tend, evdep, its_all, ite_all, Vmax, duration = event_times(outputs, dep, options)
    duration = tend - tstart

    fault_slip, evslip, rupture_length, av_slip = [],[],[],[]
    if len(tstart) > 0:
        for i in np.argsort(abs(dep)):
            z = np.sort(abs(dep))[i]
            cumslip = abs(np.array(outputs[i])[:,2])

            # -------------------- Coseismic
            Dbar = []
            for j in range(len(its_all)):
                Dbar.append(cumslip[ite_all[j]] - cumslip[its_all[j]])
                # -------------------- Event detph
                if z == evdep[j]:
                    evslip.append(cumslip[its_all[j]])
            fault_slip.append(Dbar)
        fault_slip = np.array(fault_slip)
        evslip = np.array(evslip)

        rupture_length, av_slip = analyze_events(np.sort(abs(dep)), fault_slip, tstart, options)
        
    SSE_info = {'tstart':tstart, 'tend':tend, 'its_all':its_all, 'ite_all':ite_all,\
                'evdep':evdep, 'evslip':evslip, 'fault_slip':fault_slip,\
                'duration':duration, 'Vmax':Vmax, \
                'rupture_length':rupture_length, 'av_slip':av_slip, 'options':options}
    
    return SSE_info

def event_times(outputs, dep, options):
    from scipy.signal import find_peaks
    if options['SSE'] and options['print_on']: print('Detecting SSEs')

    time = np.array(outputs[0][:,0])    
    if options['SSE'] and options['shallow']:
        di = np.where(abs(dep) <= options['Dths'])[0]
        sliprate = abs(np.array([outputs[i][:,4] for i in np.argsort(abs(dep))]))[di,:]
        z = np.sort(abs(dep))[di]
    elif options['SSE'] and options['deep']:
        di = np.where(abs(dep) >= options['Dths'])[0]
        sliprate = abs(np.array([outputs[i][:,4] for i in np.argsort(abs(dep))]))[di,:]
        z = np.sort(abs(dep))[di]
    else:
        sliprate = abs(np.array([outputs[i][:,4] for i in np.argsort(abs(dep))]))
        z = np.sort(abs(dep))

    print_string = 'Slip rate'
    unit_string = 'm/s'
    target_var = np.max(sliprate, axis=0)
    ipsr = np.argmax(sliprate, axis=0)
    if options['SSE']:
        detec_thresh = options['Vp'] * options['Vths_mult']
        thresh_seismic = options['Vths_seismic']
    else:
        detec_thresh = options['Vths']

    if options['SSE']:
        peaks,_ = find_peaks(target_var, height = (detec_thresh, thresh_seismic), distance=options['distance'], width=options['width'])
        if options['print_on']: print('Event criteria: %1.0e %s < %s < %1.0e %s'%(detec_thresh, unit_string, print_string, thresh_seismic, unit_string))
        base = np.where(target_var < detec_thresh)[0]
    else:
        peaks,_ = find_peaks(target_var, height = detec_thresh, distance=options['distance'], width=options['width'])
        if options['print_on']: print('Event criteria: %s > %1.0e %s'%(print_string, detec_thresh, unit_string))
        base = np.where(target_var < detec_thresh)[0]

    tstart, tend, evdep = [],[],[]
    no_event = False
    if len(peaks) > 0:
        its_all, ite_all = [], []
        for pks in peaks:
            append_bool = False
            if len(np.where(pks > base)[0]) > 0 and len(np.where(pks < base)[0]) == 0:
                # Possibly unfinished SSE sequence
                pot_its = base[max(np.where(pks > base)[0])]
                pot_ite = len(target_var)-1
            else:
                pot_its = base[max(np.where(pks > base)[0])]
                pot_ite = base[min(np.where(pks < base)[0])]

            # start and end seperated by more than timesteps threshold (ts_ths)
            ts_crit = pot_ite - pot_its >= options['ts_ths']
            
            # duration of the burst longer than duration threshold (dur_ths)
            dur_crit = time[pot_ite] - time[pot_its] > options['dur_ths']
            
            if ts_crit and dur_crit:
                if len(its_all) == 0:
                    append_bool = True
                elif abs(pot_its - its_all[-1]) > 0:
                    append_bool = True
                    
            if append_bool:
                its_all.append(pot_its)
                ite_all.append(pot_ite)

        if len(its_all) == 0: no_event = True
    else:
        no_event = True

    if not no_event:
        # ----- Remove events if it is only activated at specific depth: likely to be unphysical
        if options['SSE']:
            num_active_dep = np.array([np.sum(np.sum(sliprate[:, its_all[k]:ite_all[k]] > detec_thresh, axis=1) > 0) for k in range(len(its_all))])
        else:
            num_active_dep = np.array([np.sum(np.sum(sliprate[:, its_all[k]:ite_all[k]] > options['Vths'], axis=1) > 0) for k in range(len(its_all))])
        if sum(num_active_dep <= 1) > 0:
            if options['print_on']: print('Remove single-depth activated event:', np.where(num_active_dep == 1)[0])
            its_all = np.array(its_all)[num_active_dep > 1]
            ite_all = np.array(ite_all)[num_active_dep > 1]
        else:
            if options['print_on']: print('All events activate more than one depth')
        if len(its_all) == 0: no_event = True

    if not no_event and options['SSE']:
        # ----- Remove SSEs intermingles with seismic events
        peak_val_during_pick = np.array([np.max(target_var[its_all[k]:ite_all[k]]) for k in range(len(its_all))])
        if sum(peak_val_during_pick >= thresh_seismic) > 0:
            if options['print_on']: print('Remove SSEs intermingles with seismic events:', np.where(peak_val_during_pick >= thresh_seismic)[0])
            its_all = np.array(its_all)[peak_val_during_pick < thresh_seismic]
            ite_all = np.array(ite_all)[peak_val_during_pick < thresh_seismic]
        # else:
        #     if options['print_on']: print('All events far away from seismic events')
        if len(its_all) == 0: no_event = True

    if not no_event:
        evdep = np.array(z[ipsr[its_all]])
        tstart = np.array(time[its_all])
        tend = np.array(time[ite_all])
        its_all = np.array(its_all) 
        ite_all = np.array(ite_all)
        # Compute Vmax
        psr = np.max(sliprate, axis=0)
        Vmax = np.array([max(psr[its_all[i]:ite_all[i]]) for i in range(len(its_all))])
        # Compute duration
        duration = tend - tstart
        options['detec_thresh'] = detec_thresh
        if options['SSE']: options['thresh_seismic'] = thresh_seismic
    else:
        print('No seismic event detected')
        its_all, ite_all, Vmax, duration = [], [], [], []

    return tstart, tend, evdep, its_all, ite_all, Vmax, duration

# ----------

def analyze_events(dep, fault_slip, tstart, options):
    from scipy import integrate, interpolate
    if options['SSE'] and options['shallow']:
        print('Shallow SSE')
        di = np.where(abs(dep) <= options['Dths'])[0]
        z = np.sort(abs(dep))[di]
        fault_slip = fault_slip[di,:]
    elif options['SSE'] and options['deep']:
        print('Deep SSE')
        di = np.where(abs(dep) >= options['Dths'])[0]
        z = np.sort(abs(dep))[di]
        fault_slip = fault_slip[di,:]
    else:
        # print('Seismic or single SSE')
        z = np.sort(abs(dep))

    rupture_length, av_slip = [],[]    
    fault_slip = np.array(fault_slip).T         # (len(tstart), len(dep))
    if not options['SSE']: event_cluster = cluster_events(tstart)
    for ti in range(fault_slip.shape[0]):       # for each event
        fs = fault_slip[ti]
        try:
            Sths = 1e-2
            ii = np.where(fs > Sths)[0]
            if min(ii) > 0:
                ii = np.hstack(([min(ii)-1], ii))
            if max(ii) < len(fs)-1:
                ii = np.hstack((ii, [max(ii)+1]))
        except ValueError:
            Sths = 1e-3
            ii = np.where(fs > Sths)[0]
            if min(ii) > 0:
                ii = np.hstack(([min(ii)-1], ii))
            if max(ii) < len(fs)-1:
                ii = np.hstack((ii, [max(ii)+1]))
        rl = max(z[ii]) - min(z[ii])
        f = interpolate.interp1d(z, fs)
        npts = 1000
        newz = np.linspace(min(z), max(z), npts)
        Dbar = integrate.simpson(f(newz), newz)/rl
        rupture_length.append(rl)
        av_slip.append(Dbar)
    rupture_length = np.array(rupture_length)

    if not options['SSE']: 
        partial_rupture = np.where(rupture_length < options['rths'])[0]
        system_wide = np.where(rupture_length >= options['rths'])[0]

        lead_fs,major_pr,minor_pr = [],[],[]
        for k,ec in enumerate(event_cluster):
            if sum([np.logical_and(sw >= ec[0], sw <= ec[1]) for sw in system_wide]) >= 1:
                if ec[0] not in system_wide:
                    lead_fs.append(ec[0])
            elif ec[1] - ec[0] <= 4:
                minor_pr.append(k)
            else:
                major_pr.append(k)
        lead_fs = np.array(lead_fs)
        major_pr = np.array(major_pr)
        minor_pr = np.array(minor_pr)

    if not options['SSE']: 
        return rupture_length, av_slip, system_wide, partial_rupture, event_cluster, lead_fs, major_pr, minor_pr
    else: 
        return rupture_length, av_slip

def cluster_events(tstart):
    event_gap = np.diff(tstart)/yr2sec
    event_cluster = [[0,0]]
    ci = 0
    for k,eg in enumerate(event_gap):
        if eg > 1:
            event_cluster.append([k+1,k+1])
            ci += 1
        else:
            event_cluster[ci][1] = k+1
    return np.array(event_cluster)

# ----------

def compute_STF(save_dir,outputs,dep,cumslip_outputs):
    from scipy import integrate, interpolate
    tstart = cumslip_outputs['tstart']
    tend = cumslip_outputs['tend']
    params = np.load('%s/const_params.npy'%(save_dir),allow_pickle=True)
    time = np.array([outputs[i][:,0] for i in np.argsort(abs(dep))])
    sr = abs(np.array([outputs[i][:,4] for i in np.argsort(abs(dep))]))
    z = np.sort(abs(dep))*1e3
    if 'DZ' in save_dir:
        mu = params.item().get('mu_damage')*1e9
    else:
        mu = params.item().get('mu')*1e9

    npoints = 500
    f = np.array([mu * integrate.simpson(sr[:,t],z) for t in range(sr.shape[1])])
    stf = interpolate.interp1d(time[0],f)
    Fdot=np.array([stf(np.linspace(tstart[iev],tend[iev],npoints)) for iev in range(len(tstart))])
    t = np.array([np.linspace(tstart[iev],tend[iev],npoints)-tstart[iev] for iev in range(len(tstart))])
    return t,Fdot

def compute_M0(save_dir,rupture_length, av_slip,mode,Mw):
    params = np.load('%s/const_params.npy'%(save_dir),allow_pickle=True)
    if 'DZ' in save_dir:
        mu = params.item().get('mu_damage')*1e9
    else:
        mu = params.item().get('mu')*1e9
    
    rl = rupture_length * 1e3
    if mode == '1d':
        print('1D: Moment per length')
        M0 = np.array([mu * rl[iev] * av_slip[iev] for iev in range(len(av_slip))])
    elif mode == 'approx2d':
        print('Approximated 2D: Moment assuming a square fault patch')
        M0 = np.array([mu * (rl[iev]**2) * av_slip[iev] for iev in range(len(av_slip))])
    if Mw:
        print('Output in moment magnitude (Mw) instead of moment (M0)')
        return 2/3*(np.log10(M0)-9.1)
    else:
        return M0
    
def estimate_coef(x, y):
    # number of observations/points
    n = np.size(x)
 
    # mean of x and y vector
    m_x = np.mean(x)
    m_y = np.mean(y)
 
    # calculating cross-deviation and deviation about x
    SS_xy = np.sum(y*x) - n*m_y*m_x
    SS_xx = np.sum(x*x) - n*m_x*m_x
 
    # calculating regression coefficients
    b = SS_xy / SS_xx
    a = m_y - b*m_x
 
    return (a, b)

def compute_GR(save_dir,cumslip_outputs,spin_up_idx,cutoff_Mw,npts):
    rupture_length = cumslip_outputs['rupture_length']
    av_slip = cumslip_outputs['av_slip']
    Mw = compute_M0(save_dir,rupture_length, av_slip,mode='approx2d',Mw=True)
    Mw = Mw[spin_up_idx:]
    # baseval = np.linspace(min(Mw),max(Mw),npts)[:-1]
    # N = np.array([sum(Mw > mag) for mag in baseval])
    baseval = np.linspace(min(Mw),max(Mw),npts)
    # baseval = np.sort(Mw)
    N = np.array([sum(Mw >= mag) for mag in baseval])
    if cutoff_Mw == 0:
        x = baseval; y = np.log10(N)
    else:
        # ii = np.where(baseval>=cutoff_Mw)[0]; x = baseval[ii]; y = np.log10(N[ii])
        ii = np.where(baseval<=cutoff_Mw)[0]; x = baseval[ii]; y = np.log10(N[ii])
    a,b = estimate_coef(x,y)
    yN = np.power(10,a + b*x)

    return baseval,N,b,x,yN,a