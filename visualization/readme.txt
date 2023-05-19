***** Guide for using get_plots.py *****
*     By Jeena Yun (j4yun@ucsd.edu)    *
****************************************

get_plots.py is an executable python script for generating and saving plots using Tandem fault outputs.
The script requires the fault outputs to be in a single folder named 'outputs'
The script relies on individual plotting scripts: faultoutputs_vs_time.py, faultoutputs_image.py, cumslip_compute.py, cumslip_plot.py, and stress_profile.py.
See create_plot.sh for example usages.
Feel free to reach out to the author for any bugs or suggestions for improvement.

*** Required packages
numpy, matplotlib, scipy, cmcrameri (https://www.fabiocrameri.ch/colourmaps/)

*** Parameter description
[General parameters]
save_dir: Directory to output files (and to save plots). You may need to modify load/save directories.
--Wf, -wf: Maximum depth of fault. Mostly used for plotting purpose (to set ylim). [default = 24 km]
--compute, -c: If used, the raw output is processed and saved into a preferred format. Otherwise, load the processed output. Turning off after initial usage is recommended.
plot_in_sec, -sec: If used, the time axis become seconds. Otherwise, plotted in years.

[Parameters for fault output vs. time at certain depth plot]
(all optional)
--sliprate, -sr: If used, returns depth of slip rate vs. time plot. [km]
--slip, -sl: If used, returns depth of slip vs. time plot. [km]
--stress, -st: If used, returns depth of stress vs. time plot. [km]
--state_var, -sv: If used, returns depth of state_var vs. time plot. [km]

[Parameters for initial condition plot]
--stressprof, -ist: If used, returns initial stress profile used for the simulation

[Parameters for fault output image plot]
--image_sliprate, -imsr: If used, returns slip rate image plot
--image_shearT, -imst: If used, returns shear stress image plot
--image_normalT, -imnt: If used, returns normal stress image plot
--image_state_var, -imsv: If used, returns state variable image plot
--plot_in_timestep, -ts: If used, the time axis become timestep. Otherwise, plotted in time.
--vmin, -vmin: If given, min value for the colormap. Otherwise, vmin = data minimum or 1e-13 for slip rate.
--vmax, -vmax: If given, max value for the colormap. Otherwise, vmax = data maximum.

[Parameters for cumulative slip vs. depth plot]
--cumslip, -csl: If used, returns cumulative slip vs. depth returns plot.
--spin_up, -spup: If used, cumslip plot only shows slips larger than the given value (i.e., spin-up). The slip amount at given stage is removed from all following event. [m]
--dt_creep, -dtcr: Contour line interval for the creeping section. Required when cumslip plotting is on. [s]
--dt_coseismic, -dtco: Contour line interval for the coseismic section. Required when cumslip plotting is on. [s]
--dt_interm, -dtint: If given, contour line interval for the coseismic section. [s]
--Vths, -Vths: Slip-rate threshold to define coseismic section. Required when cumslip plotting is on. [m/s] [default = 1e-2 m/s]
--Vlb, -Vlb: When used with --dt_interm, becomes lower bound of slip rate of intermediate section. [m/s] [default = 1e-8 when dt_interm is non-zero; 0 otherwise]
--ev_anal, -evan: If used, returns event analyzation plot.
--depth_dist, -dd: If used, a histogram for the hypocenter depth distribution is plotted with the cumslip plot.
--rths, -rths: Rupture length threshold to define system wide event [m] [default = 10 m]
--cuttime, -ct: If given, show result up until to the given time to save computation time. Otherwise, show till the end of the simulation. [yr]
--mingap, -mg: If given, minimum seperation time between two prospective events to determine as seperate events. If the two events are seperated less than this value, they are considered as one. Used only for plotting purpose. [s] [default = 60 s]