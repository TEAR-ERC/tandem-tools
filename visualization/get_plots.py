#!/usr/bin/env python3
# By Jeena Yun (j4yun@ucsd.edu)
# Last modification: 2023.05.18.

import numpy as np
import glob
import argparse
import csv

yr2sec = 365*24*60*60
wk2sec = 7*24*60*60

# Set input parameters -------------------------------------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument("save_dir", help=": directory to output files and to save plots")
parser.add_argument("-wf","--Wf", type=float, help=": Fault length [km]")
parser.add_argument("-c","--compute", action="store_true", help=": Activate only when you want to compute")
parser.add_argument("-sec","--plot_in_sec", action="store_true", help=": Time axis in seconds",default=False)

# Fault output vs. time at certain depth
parser.add_argument("-sr","--sliprate", type=float, help=": If used, depth of slip rate vs. time plot [km]", default=0)
parser.add_argument("-sl","--slip", type=float, help=": If used, depth of slip vs. time plot [km]", default=0)
parser.add_argument("-st","--stress", type=float, help=": If used, depth of stress vs. time plot [km]", default=0)
parser.add_argument("-sv","--state_var", type=float, help=": If used, depth of state variable vs. time plot [km]", default=0)

# Input variable profile
parser.add_argument("-ist","--stressprof", action="store_true", help=": ON/OFF out stress profile")

# Fault output image
parser.add_argument("-imsr","--image_sliprate", action="store_true", help=": ON/OFF slip rate image plot")
parser.add_argument("-imst","--image_shearT", action="store_true", help=": ON/OFF shear stress image plot")
parser.add_argument("-imnt","--image_normalT", action="store_true", help=": ON/OFF normal stress image plot")
parser.add_argument("-imsv","--image_state_var", action="store_true", help=": ON/OFF state variable image plot")
parser.add_argument("-ts","--plot_in_timestep", action="store_true", help=": Time axis in timesteps",default=False)
parser.add_argument("-vmin","--vmin", type=float, help=": vmin for the plot")
parser.add_argument("-vmax","--vmax", type=float, help=": vmax for the plot")

# Cumulative slip profile related paramters
parser.add_argument("-csl","--cumslip", action="store_true", help=": ON/OFF cumulative slip profile")
parser.add_argument("-dtcr","--dt_creep", type=float, help=": Contour interval for CREEPING section [yr]")
parser.add_argument("-dtco","--dt_coseismic", type=float, help=": Contour interval for COSEISMIC section [s]")
parser.add_argument("-dtint","--dt_interm", type=float, help=": Contour interval for INTERMEDIATE section [wk]")
parser.add_argument("-Vths","--Vths", type=float, help=": Slip-rate threshold to define coseismic section [m/s]")
parser.add_argument("-Vlb","--Vlb", type=float, help=": When used with --dt_interm becomes lower bound of slip rate of intermediate section [m/s]")
parser.add_argument("-evan","--ev_anal", action="store_true", help=": ON/OFF event analyzation plot")
parser.add_argument("-dd","--depth_dist", action="store_true", help=": Plot cumslip plot with hypocenter depth distribution",default=False)
parser.add_argument("-spup","--spin_up", type=float, help=": Plot with spin-up after given slip amount",default=0)
parser.add_argument("-rths","--rths", type=float, help=": Rupture length threshold to define system wide event [m]",default=10)
parser.add_argument("-ct","--cuttime", type=float, help=": Show result up until to the given time to save computation time [yr]", default=0)
parser.add_argument("-mg","--mingap", type=float, help=": Minimum seperation time between two different events [s]", default=60)

args = parser.parse_args()

# --- Check dependencies
if args.cumslip:        # When args.cumslip are true
    if not args.dt_creep:
        parser.error('Required field \'dt_creep\' is not defined - check again')
    if not args.dt_coseismic:
        parser.error('Required field \'dt_coseismic\' is not defined - check again')
    if not args.Vths:
        print('Required field \'Vths\' not defined - using default value 1e-2 m/s')
        args.Vths = 1e-2
    dt_creep = args.dt_creep*yr2sec
    dt_coseismic = args.dt_coseismic
    if args.dt_interm:
        dt_interm = args.dt_interm*wk2sec
        if not args.Vlb:
            print('Required field \'Vlb\' not defined - using default value 1e-8 m/s')
            args.Vlb = 1e-8
    else:
        dt_interm = 0
        args.Vlb = 0

if not args.Wf:
    print('Required field Wf not defined - using default value 24 km')
    args.Wf = 24
save_dir = args.save_dir
cuttime = args.cuttime*yr2sec

# Extract data ---------------------------------------------------------------------------------------------------------------------------
if args.compute:
    print('Compute on - extract outputs...',end=' ')
    fnames = glob.glob('%s/outputs/*.csv'%(save_dir))
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
    outputs = np.load('%s/outputs.npy'%(save_dir))
    print('Load saved data: %s/outputs_depthinfo'%(save_dir))
    dep = np.load('%s/outputs_depthinfo.npy'%(save_dir))

# Fault output vs. time at certain depth -------------------------------------------------------------------------------------------------
if abs(args.sliprate)>0:
    from faultoutputs_vs_time import sliprate_time
    sliprate_time(save_dir,outputs,dep,args.sliprate,args.plot_in_sec)
    
if abs(args.slip)>0:
    from faultoutputs_vs_time import slip_time
    slip_time(save_dir,outputs,dep,args.slip,args.plot_in_sec)
    
if abs(args.stress)>0:
    from faultoutputs_vs_time import stress_time
    stress_time(save_dir,outputs,dep,args.stress,args.plot_in_sec)

if abs(args.state_var)>0:
    from faultoutputs_vs_time import state_time
    state_time(save_dir,outputs,dep,args.state_var,args.plot_in_sec)

# Fault output image ---------------------------------------------------------------------------------------------------------------------
if sum([args.image_sliprate,args.image_shearT,args.image_normalT,args.image_state_var])>0:
    from faultoutputs_image import fout_image
    if not args.vmin:
        if args.image_sliprate:
            vmin = 1e-13
        else:
            vmin = None
    else:
        vmin = args.vmin
    if not args.vmax:
        vmax = None
    else:
        vmax = args.vmax
    fout_image(args.image_sliprate,args.image_shearT,args.image_normalT,args.image_state_var,outputs,dep,args.plot_in_timestep,save_dir,args.Wf,vmin,vmax,args.plot_in_sec)

# Input variable profile -----------------------------------------------------------------------------------------------------------------
if args.stressprof:
    from stress_profile import *
    plot_stress_vs_depth(save_dir,args.Wf,outputs,dep)    
    plot_hist(save_dir,outputs,dep)

# Cumslip vs. Depth ----------------------------------------------------------------------------------------------------------------------
if args.cumslip:
    from cumslip_compute import *
    from cumslip_plot import *
    cumslip_outputs = compute_cumslip(outputs,dep,cuttime,args.Vlb,args.Vths,dt_creep,dt_coseismic,dt_interm,args.mingap)
    if args.ev_anal:
        plot_event_analyze(save_dir,args.Wf,cumslip_outputs,args.rths)

    # --- Plot the result
    if args.spin_up > 0:
        spup_cumslip_outputs = compute_spinup(outputs,dep,cuttime,cumslip_outputs,args.spin_up)
        spup_where(save_dir,args.Wf,cumslip_outputs,spup_cumslip_outputs,args.Vths,dt_coseismic,args.rths)
    else:
        spup_cumslip_outputs = None

    if args.depth_dist:
        two_set(save_dir,args.Wf,cumslip_outputs,args.Vths,dt_coseismic,args.rths,spup_cumslip_outputs)
    else:
        only_cumslip(save_dir,args.Wf,cumslip_outputs,args.Vths,dt_coseismic,args.rths,spup_cumslip_outputs)