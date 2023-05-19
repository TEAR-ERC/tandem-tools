#!/usr/bin/env python3
'''
Functions related to plotting initial stress conditions
By Jeena Yun
Last modification: 2023.05.18.
'''
import numpy as np
import matplotlib.pylab as plt

mypink = (230/255,128/255,128/255)
myblue = (118/255,177/255,230/255)
myburgundy = (214/255,0,0)
mynavy = (17/255,34/255,133/255)

# ------------------ Retrieved
def read_output(outputs,dep):
    z = np.zeros(len(dep))
    tau = np.zeros(len(dep))
    sigma = np.zeros(len(dep))

    c = 0
    for i in np.argsort(abs(dep)):
        z[c] = abs(dep[i])
        tau[c] = np.array(outputs[i])[0,3]
        sigma[c] = np.array(outputs[i])[0,5]
        c += 1
    return z,tau,sigma

# ------------------ Initial stress check
def plot_stress_vs_depth(save_dir,Wf,outputs,dep,save_on=True):
    plt.rcParams['font.size'] = '15'
    fig,ax = plt.subplots(figsize=(9,7))
    y_ret,tau,sigma = read_output(outputs,dep)

    ax.scatter(abs(tau),abs(y_ret),lw=2.5,color=myburgundy,label='Shear Stress (output)',zorder=3)
    ax.scatter(sigma,abs(y_ret),lw=2.5,color=mynavy,label='Normal Stress (output)',zorder=3)
    ax.set_xlabel('Stress [MPa]',fontsize=17)
    ax.set_ylabel('Depth [km]',fontsize=17)
    ax.set_xlim(-2,60)
    ax.set_ylim(0,Wf) # Fault depth
    ax.invert_yaxis()
    plt.grid(True)
    ax.legend(fontsize=13,loc='lower left')
    plt.tight_layout()
    if save_on:
        plt.savefig('%s/stress_profile.png'%(save_dir))

# ------------------ Histogram
def plot_hist(save_dir,outputs,dep,save_on=True):
    y_ret,tau,sigma = read_output(outputs,dep)
    plt.rcParams['font.size'] = '15'
    fig,ax=plt.subplots(ncols=2,figsize=(18,7))
    ax[0].hist(sigma,color=myblue,edgecolor='k',lw=1)
    ax[0].set_xlabel('Stress [MPa]',fontsize=17)
    ax[0].set_ylabel('Count',fontsize=17)
    ax[0].set_title('All depth',fontsize=20,fontweight='bold')

    ax[1].hist(sigma[abs(y_ret)>2],color=myblue,edgecolor='k',lw=1)
    ax[1].set_xlabel('Stress [MPa]',fontsize=17)
    ax[1].set_title('Depth > 2 km',fontsize=20,fontweight='bold')

    plt.tight_layout()
    if save_on:
        plt.savefig('%s/sigma_hist.png'%(save_dir))