#!/usr/bin/env python3
'''
Functions related to plotting cumulative slip vs. depth plot
By Jeena Yun
Last modification: 2023.07.10.
'''
import numpy as np
import matplotlib.pylab as plt
from cumslip_compute import analyze_events


mypink = (230/255,128/255,128/255)
mydarkpink = (200/255,110/255,110/255)
myblue = (118/255,177/255,230/255)
myburgundy = (214/255,0,0)
mynavy = (17/255,34/255,133/255)
mylightblue = (218/255,230/255,240/255)
myygreen = (120/255,180/255,30/255)
mylavender = (170/255,100/255,215/255)
mydarkviolet = (145/255,80/255,180/255)
pptyellow = (255/255,217/255,102/255)

yr2sec = 365*24*60*60
wk2sec = 7*24*60*60

def with_depth_dist(ax,Wf,cumslip_outputs):
    ax.hist(cumslip_outputs[1][1],bins=np.arange(0,Wf+0.2,0.2),color='k',edgecolor='k',orientation='horizontal')
    ax.set_xlabel('Counts',fontsize=30)
    ax.axes.yaxis.set_ticklabels([])
    ax.set_ylim(Wf,0)

def cumslip_basic(ax,Wf,cumslip_outputs,rths):
    # cumslip_outputs = [timeout, evout, creepout, coseisout, intermout]
    # cumslip_outputs[0] = [tstart_coseis,tend_coseis]
    # cumslip_outputs[1] = [evslip,evdep,fault_slip]
    # cumslip_outputs[2] = [cscreep,depcreep]
    # cumslip_outputs[3] = [cscoseis,depcoseis]
    # cumslip_outputs[4] = [csinterm,depinterm]
    system_wide,partial_rupture = analyze_events(cumslip_outputs,rths)[2:]

    if len(cumslip_outputs) > 4:
        ax.plot(cumslip_outputs[-1][0],cumslip_outputs[-1][1],color='yellowgreen',lw=1)
    ax.plot(cumslip_outputs[3][0],cumslip_outputs[3][1],color=mydarkpink,lw=1)
    ax.plot(cumslip_outputs[2][0],cumslip_outputs[2][1],color='0.62',lw=1)
    ax.set_ylabel('Depth [km]',fontsize=30)
    ax.set_xlabel('Cumulative Slip [m]',fontsize=30)
    if len(system_wide) > 0:
        ax.scatter(cumslip_outputs[1][0][system_wide],cumslip_outputs[1][1][system_wide],marker='*',s=700,facecolor=mydarkviolet,edgecolors='k',lw=1,zorder=3,label='System-wide events')
    if len(partial_rupture) > 0:
        ax.scatter(cumslip_outputs[1][0][partial_rupture],cumslip_outputs[1][1][partial_rupture],marker='d',s=250,facecolor=mylightblue,edgecolors='k',lw=1,zorder=3,label='Partial rupture events')
    ax.legend(fontsize=20,framealpha=1,loc='lower right')
    xl = ax.get_xlim()
    ax.set_xlim(0,xl[1])
    ax.set_ylim(0,Wf)
    ax.invert_yaxis()

def cumslip_spinup(ax,Wf,cumslip_outputs,spup_cumslip_outputs,rths):
    # spup_cumslip_outputs = [new_inits, spup_evslip, spup_cscreep, spup_cscoseis, spup_csinterm, spin_up_idx]
    # new_inits = [new_init_Sl,new_init_dp]
    system_wide,partial_rupture = analyze_events(cumslip_outputs,rths)[2:]

    if len(cumslip_outputs) > 4:
        ax.plot(spup_cumslip_outputs[4],cumslip_outputs[4][1],color='yellowgreen',lw=1)
    ax.plot(spup_cumslip_outputs[3],cumslip_outputs[3][1],color=mydarkpink,lw=1)
    ax.plot(spup_cumslip_outputs[2],cumslip_outputs[2][1],color='0.62',lw=1)
    if len(system_wide) > 0:
        ax.scatter(spup_cumslip_outputs[1][system_wide],cumslip_outputs[1][1][system_wide],marker='*',s=700,facecolor=mydarkviolet,edgecolors='k',lw=1,zorder=3,label='System-wide events')
    if len(partial_rupture) > 0:
        ax.scatter(spup_cumslip_outputs[1][partial_rupture],cumslip_outputs[1][1][partial_rupture],marker='d',s=250,facecolor=mylightblue,edgecolors='k',lw=1,zorder=3,label='Partial rupture events')
    ax.legend(fontsize=20,framealpha=1,loc='lower right')
    xl = ax.get_xlim()
    ax.set_xlim(0,xl[1])
    ax.set_ylabel('Depth [km]',fontsize=30)
    ax.set_xlabel('Cumulative Slip [m]',fontsize=30)
    ax.set_ylim(0,Wf)
    ax.invert_yaxis()

def two_set(save_dir,Wf,cumslip_outputs,Vths,dt_coseismic,rths,spup_cumslip_outputs=None,save_on=True):
    plt.rcParams['font.size'] = '27'
    fig,ax = plt.subplots(ncols=2, figsize=(24,11), gridspec_kw={'width_ratios': [4, 1]})
    plt.subplots_adjust(wspace=0.05)
    if spup_cumslip_outputs is not None:
        cumslip_spinup(ax[0],Wf,cumslip_outputs,spup_cumslip_outputs,rths)
        lab = 'spinup_'
    else:
        cumslip_basic(ax[0],Wf,cumslip_outputs,rths)
        lab = ''
    with_depth_dist(ax[1],Wf,cumslip_outputs)
    plt.tight_layout()
    if save_on:
        plt.savefig('%s/%scumslip_%d_%d_withdepth.png'%(save_dir,lab,int(Vths*100),int(dt_coseismic*10)),dpi=300)

def only_cumslip(save_dir,Wf,cumslip_outputs,Vths,dt_coseismic,rths,spup_cumslip_outputs=None,save_on=True):
    plt.rcParams['font.size'] = '27'
    fig,ax = plt.subplots(figsize=(18,11))
    if spup_cumslip_outputs is not None:
        cumslip_spinup(ax,Wf,cumslip_outputs,spup_cumslip_outputs,rths)
        lab = 'spinup_'
    else:
        cumslip_basic(ax,Wf,cumslip_outputs,rths)
        lab = ''
    plt.tight_layout()
    if save_on:
        plt.savefig('%s/%scumslip_%d_%d.png'%(save_dir,lab,int(Vths*100),int(dt_coseismic*10)),dpi=300)

def spup_where(save_dir,Wf,cumslip_outputs,spup_cumslip_outputs,Vths,dt_coseismic,rths,save_on=True):
    plt.rcParams['font.size'] = '27'
    fig,ax=plt.subplots(figsize=(19,11))
    cumslip_basic(ax,Wf,cumslip_outputs,rths)
    ax.plot(spup_cumslip_outputs[0][0],spup_cumslip_outputs[0][1],color='yellowgreen',lw=5)
    plt.tight_layout()
    if save_on:
        plt.savefig('%s/spinup_cumslip_%d_%d_where.png'%(save_dir,int(Vths*100),int(dt_coseismic*10)),dpi=300)