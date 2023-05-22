#!/usr/bin/env python3
'''
Functions related to plotting cumulative slip vs. depth plot
By Jeena Yun
Last modification: 2023.05.22.
'''
import numpy as np
import matplotlib.pylab as plt
import matplotlib.cm as pltcm
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
    xl = ax.get_xlim()
    ax.set_ylim(Wf,0)
    ax.set_xlim(xl)

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
    if len(partial_rupture) > 0:
        ax.scatter(cumslip_outputs[1][0][partial_rupture],cumslip_outputs[1][1][partial_rupture],marker='*',s=700,facecolor=mylightblue,edgecolors='k',lw=1,zorder=3,label='Partial rupture events')
    if len(system_wide) > 0:
        ax.scatter(cumslip_outputs[1][0][system_wide],cumslip_outputs[1][1][system_wide],marker='*',s=700,facecolor=mydarkviolet,edgecolors='k',lw=1,zorder=3,label='System-wide events')
    ax.legend(fontsize=25,framealpha=1,loc='lower right')
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
    if len(partial_rupture) > 0:
        ax.scatter(spup_cumslip_outputs[1][partial_rupture],cumslip_outputs[1][1][partial_rupture],marker='*',s=700,facecolor=mylightblue,edgecolors='k',lw=1,zorder=3,label='Partial rupture events')
    if len(system_wide) > 0:
        ax.scatter(spup_cumslip_outputs[1][system_wide],cumslip_outputs[1][1][system_wide],marker='*',s=700,facecolor=mydarkviolet,edgecolors='k',lw=1,zorder=3,label='System-wide events')
    ax.legend(fontsize=25,framealpha=1,loc='lower right')
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
        spin_up = True
    else:
        spin_up = False
    if not spin_up:
        cumslip_basic(ax,Wf,cumslip_outputs,rths)
    else:
        cumslip_spinup(ax,Wf,cumslip_outputs,spup_cumslip_outputs,rths)
    plt.tight_layout()
    if save_on:
        if not spin_up:
            plt.savefig('%s/cumslip_%d_%d.png'%(save_dir,int(Vths*100),int(dt_coseismic*10)),dpi=300)
        else:
            plt.savefig('%s/spinup_cumslip_%d_%d.png'%(save_dir,int(Vths*100),int(dt_coseismic*10)),dpi=300)

def spup_where(save_dir,Wf,cumslip_outputs,spup_cumslip_outputs,Vths,dt_coseismic,rths,save_on=True):
    system_wide,partial_rupture = analyze_events(cumslip_outputs,rths)[2:]
    plt.rcParams['font.size'] = '27'
    plt.figure(figsize=(19,11))

    if len(cumslip_outputs) > 4:
        plt.plot(cumslip_outputs[4][0],cumslip_outputs[4][1],color='yellowgreen',lw=1)
    plt.plot(cumslip_outputs[3][0],cumslip_outputs[3][1],color=mydarkpink,lw=1)
    plt.plot(cumslip_outputs[2][0],cumslip_outputs[2][1],color='0.62',lw=1)
    plt.plot(spup_cumslip_outputs[0][0],spup_cumslip_outputs[0][1],color='yellowgreen',lw=5)
    ev_part = plt.scatter(cumslip_outputs[1][0][partial_rupture],cumslip_outputs[1][1][partial_rupture],marker='*',s=700,facecolor=mylightblue,edgecolors='k',lw=1,zorder=3)
    ev_sys = plt.scatter(cumslip_outputs[1][0][system_wide],cumslip_outputs[1][1][system_wide],marker='*',s=700,facecolor=mydarkviolet,edgecolors='k',lw=1,zorder=3)
    plt.legend([ev_sys,ev_part],['System-wide events','Partial rupture events'],fontsize=25,framealpha=1,loc='lower right')
    plt.ylabel('Depth [km]',fontsize=30)
    plt.xlabel('Cumulative Slip [m]',fontsize=30)
    xl = plt.gca().get_xlim()
    plt.xlim(0,xl[1])
    plt.ylim(0,Wf)
    plt.gca().invert_yaxis()
    plt.tight_layout()
    if save_on:
        plt.savefig('%s/spinup_cumslip_%d_%d_where.png'%(save_dir,int(Vths*100),int(dt_coseismic*10)),dpi=300)

def plot_event_analyze(save_dir,Wf,cumslip_outputs,rths,save_on=True):
    print('Rupture length criterion:',rths,'m')
    rupture_length,av_slip,system_wide,partial_rupture = analyze_events(cumslip_outputs,rths)

    # ------ Define figure properties
    plt.rcParams['font.size'] = '15'
    plt.figure(figsize=(14,11))
    ax1 = plt.subplot2grid(shape=(2,5),loc=(0,0),colspan=3)
    plt.subplots_adjust(hspace=0.1)
    ax2 = plt.subplot2grid(shape=(2,5),loc=(1,0),colspan=3)
    ax3 = plt.subplot2grid(shape=(2,5),loc=(0,3),colspan=2,rowspan=2)
    plt.subplots_adjust(wspace=0.8)

    # cmap = pltcm.get_cmap('ocean')
    cmap = pltcm.get_cmap('gnuplot')

    # ------ Rupture length
    markers, stemlines, baseline = ax1.stem(np.arange(1,len(rupture_length)+1),rupture_length)
    plt.setp(stemlines, color='k', linewidth=2.5)
    plt.setp(markers, color='k')
    plt.setp(baseline, color='0.62')
    for i in range(len(system_wide)):
        markers, stemlines, baseline = ax1.stem(system_wide[i]+1,rupture_length[system_wide[i]])
        plt.setp(stemlines, color=cmap(np.linspace(0.5,0.9,len(system_wide))[i]), linewidth=2.6)
        plt.setp(markers, color=cmap(np.linspace(0.5,0.9,len(system_wide))[i]))
    ax1.hlines(y=rths,xmax=len(rupture_length)+1,xmin=0,color=mynavy,lw=1.5,linestyles='--')
    ax1.set_xticks(np.arange(-5,len(rupture_length)+5,5), minor=True)
    ax1.set_xlim(1-len(rupture_length)*0.05,len(rupture_length)*1.05)
    ax1.axes.xaxis.set_ticklabels([])
    ax1.set_ylabel('Rupture Length [km]',fontsize=17)
    ax1.grid(True,alpha=0.4,which='both')

    # ------ Hypocenter depth
    evdep = cumslip_outputs[1][1]
    markers, stemlines, baseline = ax2.stem(np.arange(1,len(evdep)+1),evdep)
    plt.setp(stemlines, color='k', linewidth=2.5)
    plt.setp(markers, color='k')
    plt.setp(baseline, color='0.62')
    for i in range(len(system_wide)):
        markers, stemlines, baseline = ax2.stem(system_wide[i]+1,evdep[system_wide[i]])
        plt.setp(stemlines, color=cmap(np.linspace(0.5,0.9,len(system_wide))[i]), linewidth=2.6)
        plt.setp(markers, color=cmap(np.linspace(0.5,0.9,len(system_wide))[i]))
    ax2.set_xticks(np.arange(-5,len(evdep)+5,5), minor=True)
    ax2.set_xlim(ax1.get_xlim())
    ax2.set_xlabel('Event Index',fontsize=17)
    ax2.set_ylabel('Hypocenter Depth [km]',fontsize=17)
    ax2.grid(True,alpha=0.4,which='both')

    # ------ Slip along fault for each event
    fault_z = np.array(cumslip_outputs[3][1]).T[0]
    ax3.set_prop_cycle('color',[cmap(i) for i in np.linspace(0.5,0.9,len(system_wide))])
    if len(partial_rupture) > 0:
        ax3.plot(np.array(cumslip_outputs[1][2]).T[partial_rupture[0]].T,fault_z,lw=2.5,color='0.62',label='Partial rupture events')
        if len(partial_rupture) > 1:
            ax3.plot(np.array(cumslip_outputs[1][2]).T[partial_rupture[1:]].T,np.array([fault_z for i in range(len(partial_rupture[1:]))]).T,lw=2.5,color='0.62')
    if len(system_wide) > 0:        
        ax3.plot(np.array(cumslip_outputs[1][2]).T[system_wide].T,np.array([fault_z for i in range(len(system_wide))]).T,lw=3,label=[r'Event %d ($\bar{D}$ = %2.2f m)'%(i+1,av_slip[i]) for i in system_wide])
        hyp_dep = cumslip_outputs[1][1][system_wide]
        hyp_slip = [np.array(cumslip_outputs[1][2]).T[system_wide][i][np.where(fault_z==hyp_dep[i])[0][0]] for i in range(len(system_wide))]
        ax3.scatter(hyp_slip,hyp_dep,marker='*',s=300,facecolor=mydarkviolet,edgecolors='k',lw=1,zorder=3,label='Hypocenter')
    if len(partial_rupture) > 0 or len(system_wide) > 0:
        ax3.legend(fontsize=13)
    xl = ax3.get_xlim()
    ax3.set_xlabel('Slip [m]',fontsize=17)
    ax3.set_ylabel('Depth [km]',fontsize=17)
    ax3.set_xlim(xl)
    ax3.set_ylim(0,Wf)
    ax3.invert_yaxis()
    ax3.grid(True,alpha=0.4)

    plt.tight_layout()
    if save_on:
        plt.savefig('%s/analyze_events.png'%(save_dir),dpi=300)