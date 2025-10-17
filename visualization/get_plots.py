#!/usr/bin/env python3
'''
An executable plotting script for Tandem to generate figures from Tandem output
Author: Jeena Yun (j4yun@ucsd.edu)
Last modification: 2025.10.15.
'''
# ------------------------- PLEASE CHANGE THIS LINE BEFORE USING ------------------------- #
root_dir = '/ROOT/DIRECTORY/TO/YOUR/TANDEM/OUTPUTS'
# ---------------------------------------------------------------------------------------- #

import argparse
from read_outputs import OUTPUTS

yr2sec = 365*24*60*60
wk2sec = 7*24*60*60


# ----- Set input parameters
parser = argparse.ArgumentParser()
parser.add_argument("job_name", help=": Name of the Tandem job")

# Output prefix
parser.add_argument("-fp_prefix", "--faultp_prefix", help=": Prefix of fault probe outputs")
parser.add_argument("-dp_prefix", "--domainp_prefix", help=": Prefix of domain probe outputs")

# Plot type
parser.add_argument("-im", "--image", type=str, choices=['state','slip','sliprate','shearT','delshearT','normalT','delnormalT','dCFS'], help=": Type of image plot ['state','slip','sliprate','shearT','delshearT','normalT','delnormalT']")
parser.add_argument("-csl", "--cumslip", action="store_true", help=": ON/OFF cumulative slip profile")
parser.add_argument("-ts", "--timeseries", type=str, choices=['state','slip','sliprate','stress','displacement'], help=": Plot time series of the given variable")
parser.add_argument("-loc", "--target_loc", nargs='+', type=float, help=": When used with --timeseries, plot timeseries at the given receiver location. For displacement, it takes (x,y) location, while on-fault variables only need depth. If not given, plot peak along the fault.")

# Time control
parser.add_argument("-step", "--plot_in_timestep", action="store_true", help=": Time axis in timesteps", default=False)
parser.add_argument("-sec", "--plot_in_sec", action="store_true", help=": Time axis in seconds",default=False)

args = parser.parse_args()

kwargs = {k: v for k, v in vars(args).items() if v is not None}

# ----- Define output directory
save_dir = root_dir + args.job_name
kwargs['save_dir'] = save_dir

# ------ Extract output
var = OUTPUTS(save_dir)
var.get_outputs(abs_on=True, faultp_prefix=args.faultp_prefix)
outputs, dep = var.outputs, var.dep

# ------ Compute event info
if (args.cumslip or args.image):
    event_info = var.get_event_info(**kwargs, save_on=True)
else:
    event_info = []

# ------ Fault output image 
if args.image:
    from faultoutputs_image import fout_image    
    fout_image(outputs, dep, event_info, **kwargs)

# ----- Cumslip vs. Depth
if args.cumslip:
    from cumslip_plot import only_cumslip
    from cumslip_compute import compute_cumslip
    
    cumslip_outputs = compute_cumslip(outputs, dep, event_info, **kwargs)
    # --- Plot the result
    only_cumslip(cumslip_outputs, event_info, **kwargs)

# ----- Timeseries
if args.timeseries:    
    if args.timeseries == 'displacement' or args.timeseries == 'stress':
        import matplotlib.pylab as plt
        plt.rcParams['font.size'] = '11'
        fig, [ax,ax2] = plt.subplots(nrows=2, figsize=(10,9), sharex=True)
        if args.timeseries == 'displacement':
            if not args.target_loc:
                raise SyntaxError('No location provided - please provide location for displacement plot')
            if len(args.target_loc) == 1:
                raise SyntaxError('Only x or y coordinate received - please provide a (x, y) location for displacement plot')
            var.get_domain_outputs(domainp_prefix=args.domainp_prefix)
            var.timeseries(ax, 'u0', **kwargs)
            ax.set_xlabel('')
            var.timeseries(ax2, 'u1', target_loc=args.target_loc, plot_in_sec=args.plot_in_sec)
            figname = 'displacement_at_%dm_%dm'%(args.target_loc[0]*1e3, args.target_loc[1]*1e3)
            fig_title = 'At (%1.1f km, %1.1f km)'%(args.target_loc[0], args.target_loc[1])
        if args.timeseries == 'stress':
            if not args.target_loc:
                kwargs['target_loc'] = None
                figname = 'peak_stresses'
            elif len(args.target_loc) == 2:
                raise SyntaxError('Two values received - please provide only depth for on-fault variables plot')
            else:
                figname = 'stresses_at_%dm_depth'%(abs(args.target_loc[0])*1e3)
                if args.target_loc[0] < 1e-1:
                    fig_title = 'At surface'
                else:
                    fig_title = 'At %1.1f km Depth'%(abs(args.target_loc[0]))
            var.timeseries(ax, 'shearT', target_loc=args.target_loc, plot_in_sec=args.plot_in_sec)
            ax.set_xlabel('')
            var.timeseries(ax2, 'normalT', target_loc=args.target_loc, plot_in_sec=args.plot_in_sec)
        if args.target_loc:
            fig.suptitle(fig_title,fontsize=15,fontweight = 'bold')
        plt.tight_layout()
        if args.plot_in_timestep: plt.savefig('%s/%s_timesteps.png'%(save_dir, figname), dpi=150)
        else: plt.savefig('%s/%s.png'%(save_dir, figname), dpi=150)
    else:
        if not args.target_loc:
            kwargs['target_loc'] = None
        elif len(args.target_loc) == 2:
            raise SyntaxError('Two values received - please provide only depth for on-fault variables plot')
        var.timeseries(None, args.timeseries, **kwargs)
    