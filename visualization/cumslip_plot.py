#!/usr/bin/env python3
'''
Functions related to plotting cumulative slip vs. depth plot
By Jeena Yun
Last modification: 2025.10.17.
'''
import numpy as np
import matplotlib.pylab as plt


mydarkpink = (200/255,110/255,110/255)
mydarkviolet = (145/255,80/255,180/255)
yr2sec = 365*24*60*60
wk2sec = 7*24*60*60

def only_cumslip(cumslip_outputs, event_info, **kwargs):
    options = {
        'dt_creep' : 5*yr2sec,
        'dt_coseismic' : 1,
        'print_on' : True
    }
    options.update(kwargs)

    Wf = max([np.max(abs(np.array(cumslip_outputs['depcreep']))), np.max(abs(np.array(cumslip_outputs['depcoseis'])))])

    plt.rcParams['font.size'] = '27'
    fig,ax = plt.subplots(figsize=(18,11))

    ax.plot(cumslip_outputs['cscoseis'],cumslip_outputs['depcoseis'],color=mydarkpink,lw=1)
    ax.plot(cumslip_outputs['cscreep'],cumslip_outputs['depcreep'],color='0.62',lw=1)
    ax.set_ylabel('Depth [km]',fontsize=30)
    ax.set_xlabel('Cumulative Slip [m]',fontsize=30)
    ax.scatter(event_info['evslip'],event_info['evdep'],marker='*',s=700,facecolor=mydarkviolet,edgecolors='k',lw=1,zorder=3,label='Earthquakes')
    ax.legend(fontsize=20,framealpha=1,loc='lower right')
    xl = ax.get_xlim()
    ax.set_xlim(0, xl[1])
    ax.set_ylim(Wf, 0)
    plt.tight_layout()
    
    plt.savefig('%s/cumslip.png'%(options['save_dir']),dpi=300)