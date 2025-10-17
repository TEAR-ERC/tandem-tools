#!/usr/bin/env python3
'''
Functions related to reading tandem outputs: fault/fault_probe/domain_probe
By Jeena Yun
Last modification: 2025.03.14.
'''
import numpy as np
from glob import glob
import os

yr2sec = 365*24*60*60

def base_event_criteria():
    options = {
        'mode' : 'sliprate',
        'Pths_mult' : 1e-2,
        'Vths' : 1e-2,
        'distance' : 500,
        'width' : 10,
        'ts_ths' : 5,
        'dur_ths' : 1,
        'rths' : 10,
        'SSE' : False,
        'print_on' : True,
        'save_on' : True }
    return options

def read_fault_probe_outputs(save_dir, options):
    import time
    if os.path.exists('%s/outputs.npy'%(save_dir)):
        # outputs, dep = load_fault_probe_outputs(save_dir, options)
        if options['print_on']: print('Load saved data: %s/outputs.npy'%(save_dir))
        outputs = np.load('%s/outputs.npy'%(save_dir))
        if options['print_on']: print('Load saved data: %s/outputs_depthinfo.npy'%(save_dir))
        dep = np.load('%s/outputs_depthinfo.npy'%(save_dir))
    else:
        import pandas as pd
        # --- Read outputs
        fnames = glob('%s/%s*.csv'%(save_dir, options['faultp_prefix']))
        if len(fnames) == 0:
            raise NameError('No fault probe output found in %s - check the input'%(save_dir))
        
        # --- Start computing output
        if options['print_on']: print('Start computing output... ',end='')
        ti = time.time()
        len_dat = 0
        outputs,dep=[],[]
        for fname in np.sort(fnames):
            # dat = pd.read_csv(fname,delimiter=',',skiprows=1,nrows=200216) # only for trial_11
            dat = pd.read_csv(fname,delimiter=',',skiprows=1)
            # ** Sanitary check 1: check if the each csv file has the same size
            if len_dat == 0:
                len_dat = dat.values.shape[0]
            elif dat.values.shape[0] != len_dat:
                raise ValueError('The size of the output is not uniform: %d -> %d'%(len_dat,dat.values.shape[0]))
            stloc = pd.read_csv(fname,nrows=1,header=None)
            dep.append(float(stloc.values[0][-1].split(']')[0]))
            outputs.append(dat.values)
            len_dat = dat.values.shape[0]
        outputs = np.array(outputs)
        dep = np.array(dep)
        # ** Sanitary check 2: check if there's any non-numeric values
        if np.isnan(outputs).any(): raise ValueError('At least one non-numeric values detected: check your output')
        if options['print_on']: print('Done! (%2.4f s)'%(time.time()-ti))

        # --- Sort by depth
        ii = np.argsort(abs(dep))
        outputs = outputs[ii]
        dep = dep[ii]
        if options['save_on']:
            print('Save data...',end=' ')
            np.save('%s/outputs'%(save_dir),outputs)
            np.save('%s/outputs_depthinfo'%(save_dir),dep)
            print('done!')
    return outputs,dep

def read_domain_probe_outputs(save_dir, options):
    if os.path.exists('%s/domain_probe_outputs.npy'%(save_dir)):
        # outputs,xyloc = load_domain_probe_outputs(save_dir)
        if options['print_on']: print('Load saved data: %s/domain_probe_outputs.npy'%(save_dir))
        outputs = np.load('%s/domain_probe_outputs.npy'%(save_dir))
        if options['print_on']: print('Load saved data: %s/domain_probe_outputs_xyloc.npy'%(save_dir))
        xyloc = np.load('%s/domain_probe_outputs_xyloc.npy'%(save_dir))
    else:
        import pandas as pd
        fnames = glob('%s/%s*.csv'%(save_dir, options['domainp_prefix']))
        if len(fnames) == 0:
            raise NameError('No domain probe output found in %s - check the input'%(save_dir))
        outputs,xyloc=[],[]
        for fname in np.sort(fnames):
            dat = pd.read_csv(fname,delimiter=',',skiprows=2)
            stloc = pd.read_csv(fname,nrows=1,header=None)
            xyloc.append([float(stloc.values[0][0].split('[')[-1]),float(stloc.values[0][-1].split(']')[0])])
            outputs.append(dat.values)
        outputs = np.array(outputs)
        xyloc = np.array(xyloc)
        if options['save_on']:
            print('Save data...',end=' ')
            np.save('%s/domain_probe_outputs'%(save_dir),outputs)
            np.save('%s/domain_probe_outputs_xyloc'%(save_dir),xyloc)
            print('done!')
    return outputs,xyloc

def load_checkpoint_info(save_dir):
    import pandas as pd
    dat = np.sort(pd.read_csv('%s/checkpoint_info.csv'%(save_dir),delimiter=',').values, axis=0)
    return dat

def load_short_fault_probe_outputs(save_dir, indx, print_on=True):
    if print_on: print('Load saved data: %s/short_outputs_%d'%(save_dir,indx))
    outputs = np.load('%s/short_outputs_%d.npy'%(save_dir,indx))
    if print_on: print('Load saved data: %s/outputs_depthinfo'%(save_dir))
    dep = np.load('%s/outputs_depthinfo.npy'%(save_dir))
    return outputs,dep

def get_event_info(outputs, dep, **kwargs):
    options = base_event_criteria()
    options.update(kwargs)
    if options['print_on']: print('rths = %d km'%(options['rths']))

    branch_name = options['save_dir'].split('/')[-1]
    cumslip_dir = options['save_dir']
    if options['mode'] == 'potency': decor = 'Pths_mult'
    elif options['mode'] == 'sliprate': decor = 'Vths'

    fname = '%s/event_info_%s_%1.0e_rths_%d_width_%d_dist_%d.npy'%(cumslip_dir, decor, options['Pths_mult'], options['rths'], options['width'], options['distance'])
    exist_cumslip = os.path.exists(fname)
    if not exist_cumslip or 'pert' in branch_name:
        from event_analyze import detect_events
        print('Compute event info')
        event_info = detect_events(outputs, dep, options)
    else:
        if options['print_on']: print('Load saved data: %s'%(fname))
        event_info = np.load('%s'%(fname), allow_pickle=True).item()
    return event_info


def read_pvd(fname):
    if not os.path.exists(fname):
        raise NameError('No file %s found - check the input'%(fname))
    fid = open(fname,'r')
    lines = fid.readlines()
    time = []
    for line in lines:
        if line.split('<')[1].split()[0] == 'DataSet':
            time.append(float(line.split('\"')[1]))
    fid.close()
    time = np.array(time)
    return time

def fault_output_idx_lab(target_var):
    idx = {'time': 0, 'state': 1, 'slip': 2, 'shearT': 3, 'sliprate': 4, 'normalT': 5}
    label = {'time': 'Time [s]', \
        'state': r'State Variable $\psi$', \
        'slip': 'Slip [m]', \
        'shearT': 'Shear Traction [MPa]', \
        'sliprate': 'Slip Rate [m/s]', \
        'normalT': 'Normal Stress [MPa]'} 
    return idx[target_var], label[target_var]

def domain_output_idx_lab(target_var):
    idx = {'time': 0, 'u0': 1, 'u1': 2}
    label = {'time': 'Time [s]', \
        'u0': 'Horizontal Displacement [m]', \
        'u1': 'Vertical [m]'} 
    return idx[target_var], label[target_var]
    
def read_csv(fname,target_var):
    import pandas as pd
    outputs = pd.read_csv(fname,delimiter=',',skiprows=1).values
    stloc = pd.read_csv(fname,nrows=1,header=None)
    dep = float(stloc.values[0][-1].split(']')[0])
    time = outputs[:,0]
    idx,lab = fault_output_idx_lab(target_var)
    var = outputs[:,idx]
    return dep, time, var, lab

def min_dist(_x,_y,d,z):
    min_dist = np.min(np.sqrt(np.square(d-_x) + np.square(z-_y)))
    return min_dist

# ----------------
class OUTPUTS:
    def __init__(self, save_dir):
        self.yr2sec = yr2sec
        self.save_dir = save_dir

    def get_outputs(self, **kwargs):
        options = {
            'save_on' : False,
            'abs_on' : False,
            'ckp_load' : False,
            'theta' : False,
            'print_on' : True
        }
        options.update(kwargs)
        self.outputs, self.dep = read_fault_probe_outputs(self.save_dir, options)
        if options['ckp_load']: self.outputs = self.outputs[:,1:,:]
        for target_var in ['time', 'state', 'slip', 'shearT', 'sliprate', 'normalT']:
            idx,_ = fault_output_idx_lab(target_var)
            if target_var == 'time':
                self.time = self.outputs[0,:,idx]
            elif options['abs_on'] and target_var in ['shearT', 'sliprate']: 
                setattr(self, target_var, abs(self.outputs[:,:,idx]))
            elif target_var == 'state':
                self.state = self.outputs[:,:,idx]
            else:
                setattr(self, target_var, self.outputs[:,:,idx])

    def get_domain_outputs(self, **kwargs):
        options = {
            'save_on' : False,
            'print_on' : True
        }
        options.update(kwargs)
        self.dp_outputs, self.xyloc = read_domain_probe_outputs(self.save_dir, options)
        for target_var in ['time', 'u0', 'u1']:
            idx,_ = domain_output_idx_lab(target_var)
            if target_var == 'time':
                self.time = self.outputs[0,:,idx]
            else:
                setattr(self, target_var, self.outputs[:,:,idx])

    def get_event_info(self, **kwargs):
        if hasattr(self, 'outputs'):
            self.event_info = get_event_info(self.outputs, self.dep, **kwargs)
        else:
            self.event_info = get_event_info([], [], **kwargs)
        return self.event_info

    def read_fp_csv(self, csvname, target_var):
        fname = '%s/%s'%(self.save_dir,csvname)
        self.fp_dep, self.fp_time, self.fp_data, self.fp_lab = read_csv(fname,target_var)

    # ----- Plots
    def timeseries(self, ax, target_var, **kwargs):
        import myplots
        mp = myplots.Figpref()
        options = mp.default_options.copy()
        options.update(kwargs)
        create_fig = False
        peak_mode = False

        if target_var in ['u0', 'u1']:
            loc_info = self.xyloc
        else:
            loc_info = self.dep

        if ax is None:
            import matplotlib.pylab as plt
            plt.rcParams['font.size'] = '11'
            fig,ax = plt.subplots(figsize=(10,5))
            create_fig = True
        
        # --- Set x-axis scale and x-label
        time, options['xlab'] = mp.set_time(self.time, options['plot_in_sec'])

        # --- 
        if options['target_loc'] is not None:
            if len(options['target_loc']) == 1:
                indx = np.argmin(abs(abs(loc_info) - abs(options['target_loc'][0])))
                if 'T' in target_var:
                    options['title_on'] = False
                else:
                    if options['target_loc'][0] < 1e-1:
                        options['fig_title'] = 'At surface'
                        figname = '%s_at_surface'%(target_var)
                    else:
                        options['fig_title'] = 'At %1.1f km Depth'%abs(self.dep[indx])
                        figname = '%s_at_%dm_depth'%(target_var, abs(options['target_loc'][0])*1e3)
                    if options['print_on']: print(options['fig_title'])
            elif len(options['target_loc']) == 2:
                dist = np.square(loc_info[0] - options['target_loc'][0]) + np.square(loc_info[1] - options['target_loc'][1])
                indx = np.argmin(dist)
                options['title_on'] = False
            var = getattr(self, target_var)[indx, :]
        else:
            if options['print_on']: print('Peak Mode')
            peak_mode = True
            figname = 'peak_%s'%(target_var)
            var = np.max(abs(getattr(self, target_var)), axis=0)

        # --- Set y-scale and y-label
        _, ylab = fault_output_idx_lab(target_var)
        if peak_mode:
            ylab = 'Peak %s'%(ylab)
        if options['log10'] and target_var == 'sliprate':
            if options['print_on']: print('Slip rate - log10 scale applied')
            ylab = '$\log_{10}$(%s)'%(ylab)
            var = np.log10(var)
        options['ylab'] = ylab
        
        # --- Plot the time series
        ax = mp.plot_timeserise(ax, time, var, **options)
        if create_fig: 
            plt.tight_layout()
            if options['save_on']:
                # if peak_mode:
                #     figname = 'peak_%s'%(target_var)
                # elif len(options['target_loc']) == 1:
                #     figname = '%s_at_%dm_depth'%(target_var, abs(options['target_loc'][0])*1e3)
                # elif len(options['target_loc']) == 2:
                #     figname = '%s_at_%dm_%dm'%(target_var, abs(options['target_loc'][0])*1e3, abs(options['target_loc'][1])*1e3)
                if options['plot_in_timestep']: plt.savefig('%s/%s_timesteps.png'%(options['save_dir'], figname), dpi=150)
                else: plt.savefig('%s/%s.png'%(options['save_dir'], figname), dpi=150)
        return ax