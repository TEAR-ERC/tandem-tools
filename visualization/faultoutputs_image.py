#!/usr/bin/env python3
'''
Functions related to plotting spatio-temporal evolution of variables as an image
By Jeena Yun
Last modification: 2023.07.10.
'''
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import cmcrameri.cm as cram
import warnings

yr2sec = 365*24*60*60
wk2sec = 7*24*60*60

mypink = (230/255,128/255,128/255)
myblue = (118/255,177/255,230/255)
myburgundy = (214/255,0,0)
mynavy = (17/255,34/255,133/255)


# fields: Time [s] | state [s?] | cumslip0 [m] | traction0 [MPa] | slip-rate0 [m/s] | normal-stress [MPa]
# Index:     0     |      1     |       2      |        3        |         4        |          5 

def fout_image(lab,outputs,dep,cumslip_outputs,save_dir,Wf,
               rths,vmin,vmax,Vths,zoom_frame,plot_in_timestep=True,plot_in_sec=False,save_on=True):
    from cumslip_compute import analyze_events
    system_wide,partial_rupture = analyze_events(cumslip_outputs,rths)[2:]
    processed_outputs = class_figtype(zoom_frame,outputs,cumslip_outputs)[0]
    X,Y,var = get_var(lab,processed_outputs,dep,plot_in_timestep,plot_in_sec)
    if vmin is None:
        vmin = np.min(var)
    if vmax is None:
        vmax = np.max(var)
    plot_image(X,Y,var,lab,outputs,cumslip_outputs,system_wide,partial_rupture,save_dir,Wf,zoom_frame,plot_in_timestep,vmin,vmax,Vths,plot_in_sec,save_on)

def get_var(lab,outputs,dep,plot_in_timestep,plot_in_sec):
    if lab == 'sliprate':
        idx = 4
    elif lab == 'shearT':
        idx = 3
    elif lab == 'normalT':
        idx = 5
    elif lab == 'state_var':
        idx = 1
    var = np.array([outputs[i][:,idx] for i in np.argsort(abs(dep))])
    if lab == 'shearT' or lab == 'normalT':
        var = abs(var)

    if plot_in_timestep:
        print('Plot in time steps')
        xax = np.arange(var.shape[1])
    elif plot_in_sec:
        print('Plot in time [s]')
        xax = np.array(outputs[0][:,0])
    else:
        print('Plot in time [yrs]')
        xax = np.array(outputs[0][:,0]/yr2sec)
    X,Y = np.meshgrid(xax,np.sort(abs(dep)))

    return X,Y,var

def gen_cmap(lab,vmin,vmax,Vths):
    if lab == 'sliprate':
        cb_label = 'Slip Rate [m/s]'
        cm = mpl.colormaps['RdYlBu_r']
        col_list = [cm(i)[0:3] for i in [0.15,0.5,0.8,0.9]]
        col_list = [cm(0.15)[0:3],mpl.colormaps['jet'](0.67),cm(0.8)[0:3],cm(0.9)[0:3]]
        col_list.insert(0,(0,0,0))
        col_list.insert(3,mpl.colors.to_rgb('w'))
        col_list.append(mpl.colormaps['turbo'](0))        
        V0 = 1e-6
        Vp = 1e-9
        # V0 = float(input('Reference slip velocity V0 [m/s]: '))
        # Vp = float(input('Plate velocity Vp [m/s]: '))
        float_list = [0,mpl.colors.LogNorm(vmin,vmax)(Vp),mpl.colors.LogNorm(vmin,vmax)(V0)]
        [float_list.append(k) for k in np.linspace(mpl.colors.LogNorm(vmin,vmax)(Vths),1,4)]
        cmap_n = get_continuous_cmap(col_list,input_hex=False,float_list=float_list)
    elif lab == 'shearT':
        cb_label = 'Shear Traction Change [MPa]'
        cm = cram.vik
        col_list = [cm(i) for i in np.linspace(0,1,6)]
        col_list.insert(3,mpl.colors.to_rgb('w'))
        float_list = [mpl.colors.Normalize(vmin,vmax)(i) for i in np.linspace(vmin,0,4)]
        [float_list.append(mpl.colors.Normalize(vmin,vmax)(k)) for k in np.linspace(0,vmax,4)[1:]]
        cmap_n = get_continuous_cmap(col_list,input_hex=False,float_list=float_list)
    elif lab == 'normalT':
        cb_label = 'Normal Stress [MPa]'
        cmap_n = cram.davos
    elif lab == 'statevar':
        cb_label = 'State Variable [1/s]'
        cmap_n = 'magma'
    return cmap_n,cb_label

def class_figtype(zoom_frame,outputs,cumslip_outputs,print_on=True):
    tstart,tend,evdep = cumslip_outputs[0][0],cumslip_outputs[0][1],cumslip_outputs[1][1]
    its_all = np.array([np.argmin(abs(outputs[0][:,0]-t)) for t in tstart])
    ite_all = np.array([np.argmin(abs(outputs[0][:,0]-t)) for t in tend])
    if len(zoom_frame) == 0: # Full image
        if print_on: print('Full image')
        processed_outputs = outputs
        tsmin,tsmax,its,ite,buffer1,buffer2,iev1,iev2 = [],[],[],[],[],[],[],[]
        xl_opt = 1
        scatter_opt = 1
        vlines_opt = 0
        txt_opt = 0
        xlab_opt = 1
        name_opt = 1
    elif len(zoom_frame) == 2: # Zoom in of full image
        if print_on: print('Zoom in of full image')
        tsmin,tsmax = zoom_frame[0],zoom_frame[1]
        if tsmin > outputs.shape[1]:
            ValueError('tmin > max. timestep - check input')
        elif tsmax > outputs.shape[1]:
            warnings.warn('tsmax > max. timestep')
        its_all = np.array([np.argmin(abs(outputs[0][:,0]-t)) for t in tstart]) - tsmin
        processed_outputs = outputs[:,tsmin:tsmax,:]
        its,ite,buffer1,buffer2,iev1,iev2 = [],[],[],[],[],[]
        xl_opt = 1
        scatter_opt = 1
        vlines_opt = 0
        txt_opt = 1
        xlab_opt = 2
        name_opt = 2
    else:
        buffer1,buffer2 = zoom_frame[-2],zoom_frame[-1]
        xl_opt = 2
        scatter_opt = 2
        xlab_opt = 1
        if len(zoom_frame) == 3: # Single coseismic event
            if print_on: print('Single coseismic event')
            iev1,iev2 = zoom_frame[0],[]
            its,ite = its_all[iev1],ite_all[iev1]
            vlines_opt = 1
            txt_opt = 2
            name_opt = 3
        elif len(zoom_frame) == 4 and zoom_frame[0]>=0: # Multiple coseismic event
            if print_on: print('Multiple coseismic events')
            iev1,iev2 = zoom_frame[0],zoom_frame[1]
            its,ite = its_all[iev1],ite_all[iev2]
            vlines_opt = 2
            txt_opt = 4
            name_opt = 4
        elif len(zoom_frame) == 4 and zoom_frame[0]<0: # Interseismic event
            if print_on: print('Interseismic period')
            iev1,iev2 = abs(zoom_frame[0]),abs(zoom_frame[1])
            its,ite = ite_all[iev1],its_all[iev2]
            vlines_opt = 1
            txt_opt = 3
            name_opt = 5
        processed_outputs = outputs[:,its-buffer1:ite+buffer2,:]
        tsmin,tsmax = [],[]
    fig_opts = [xl_opt,scatter_opt,vlines_opt,txt_opt,xlab_opt,name_opt]
    return processed_outputs,evdep,its_all,ite_all,tsmin,tsmax,its,ite,buffer1,buffer2,iev1,iev2,fig_opts

def decoration(time,zoom_frame,outputs,cumslip_outputs,Wf,acolor,system_wide,partial_rupture,plot_in_timestep,plot_in_sec):
    evdep,its_all,ite_all,tsmin,tsmax,its,ite,buffer1,buffer2,iev1,iev2,fig_opts = class_figtype(zoom_frame,outputs,cumslip_outputs,print_on=False)[1:]    
    xl_opt,scatter_opt,vlines_opt,txt_opt,xlab_opt,name_opt = fig_opts

    if xl_opt == 1:
        xl = plt.gca().get_xlim()
    elif xl_opt == 2:
        width = np.round(time[ite]-time[its]) # or width = np.round(X[0][-buffer2]-X[0][buffer1])
        xl = [time[its]-width/6,time[ite]+width/6]

    if scatter_opt == 1:
        if len(system_wide) > 0:
            plt.scatter(its_all[system_wide],evdep[system_wide],s=200,marker='*',facecolor='w',edgecolor='k',lw=1.5,zorder=3,label='Full rupture events')
        if len(partial_rupture) > 0:
            plt.scatter(its_all[partial_rupture],evdep[partial_rupture],s=100,marker='d',facecolor='w',edgecolor='k',lw=1.5,zorder=3,label='Partial rupture events')
    elif scatter_opt == 2:
        if len(system_wide) > 0:
            plt.scatter(time[its_all][system_wide],evdep[system_wide],s=200,marker='*',facecolor='w',edgecolor='k',lw=1.5,zorder=3,label='Full rupture events')
        if len(partial_rupture) > 0:
            plt.scatter(time[its_all][partial_rupture],evdep[partial_rupture],s=100,marker='d',facecolor='w',edgecolor='k',lw=1.5,zorder=3,label='Partial rupture events')
    plt.legend(fontsize=20,framealpha=1,loc='lower right')

    if vlines_opt == 1:
        plt.vlines(x=[time[its],time[ite]],ymin=0,ymax=Wf,linestyles='--',color=acolor,lw=1.5)
    elif vlines_opt == 2:
        plt.vlines(x=time[its_all[iev1:iev2+1]],ymin=0,ymax=Wf,linestyles='--',color=acolor,lw=1.5)

    if txt_opt == 1:
        for k in range(len(its_all)):
            if its_all[k]>=0 and its_all[k]<=tsmax-tsmin:
                plt.text(its_all[k],evdep[k]-0.2,'%d'%(k),color=acolor,fontsize=20,ha='right',va='bottom')
    elif txt_opt == 2:
        plt.text(time[its]+width/18,23,'Event %d'%(iev1),fontsize=30,fontweight='bold',color=acolor,ha='left',va='bottom') # coseismic
    elif txt_opt == 3:
        plt.text(time[its]-width/150,21.5,'Event %d'%(iev1),fontsize=30,fontweight='bold',color=acolor,ha='right',va='bottom',rotation=90)
        plt.text(time[ite]+width/150,21.5,'Event %d'%(iev2),fontsize=30,fontweight='bold',color=acolor,ha='left',va='bottom',rotation=90)
    elif txt_opt == 4:
        for k,evts in enumerate(its_all):
            if evts>=its-buffer1 and evts<=ite+buffer2:
                plt.text(time[evts],evdep[k]-0.2,'%d'%(k),color=acolor,fontsize=20,ha='right',va='bottom')

    if xl_opt == 1:
        plt.xlim(0,xl[1])
    elif xl_opt == 2:
        plt.xlim(xl)
        
    plt.ylim(0,Wf)
    plt.gca().invert_yaxis()
    plt.ylabel('Depth [km]',fontsize=30)

    if xlab_opt == 1:
        if plot_in_timestep:
            plt.xlabel('Timesteps',fontsize=30)
        elif plot_in_sec:        
            plt.xlabel('Time [s]',fontsize=30)
        else:        
            plt.xlabel('Time [yrs]',fontsize=30)
    elif xlab_opt == 2:
        if plot_in_timestep:
            plt.xlabel('Timesteps - %d'%(tsmin),fontsize=30)
        elif plot_in_sec:        
            plt.xlabel('Time [s] - %2.2f'%(time[tsmin]),fontsize=30)
        else:        
            plt.xlabel('Time [yrs] - %2.2f'%(time[tsmin]/yr2sec),fontsize=30)

    if name_opt == 1:
        fig_name = '_image'
    elif name_opt == 2:
        num1,num2 = tsmin,tsmax
        while num2 >= 10:
            num1 /= 10
            num2 /= 10
        fig_name = '_zoom_image_%d_%d'%(num1,num2)
    elif name_opt == 3:
        fig_name = '_zoom_image_ev%d'%(iev1)
    elif name_opt == 4:
        fig_name = '_zoom_image_ev%dto%d'%(iev1,iev2)
    elif name_opt == 5:
        fig_name = '_zoom_image_interevent_%d_%d'%(iev1,iev2)

    return fig_name

def plot_image(X,Y,var,lab,outputs,cumslip_outputs,system_wide,partial_rupture,save_dir,Wf,zoom_frame,plot_in_timestep,vmin,vmax,Vths,plot_in_sec,save_on):
    plt.rcParams['font.size'] = '27'
    plt.figure(figsize=(20.6,11))

    cmap_n,cb_label = gen_cmap(lab,vmin,vmax,Vths)

    if lab == 'sliprate':
        cb = plt.pcolormesh(X,Y,var,cmap=cmap_n,norm=mpl.colors.LogNorm(vmin,vmax))
        acolor = 'w'
    elif lab == 'shearT':
        var = np.array([var[i,:]-var[i,0] for i in range(var.shape[0])])
        cb = plt.pcolormesh(X,Y,var,cmap=cmap_n,vmin=vmin,vmax=vmax)
        acolor = 'k'
    else:
        cb = plt.pcolormesh(X,Y,var,cmap=cmap_n,vmin=vmin,vmax=vmax)
        acolor = 'w'        
    plt.colorbar(cb,extend='both').set_label(cb_label,fontsize=30,rotation=270,labelpad=30)

    time = outputs[0,:,0]
    fig_name = decoration(time,zoom_frame,outputs,cumslip_outputs,Wf,acolor,system_wide,partial_rupture,plot_in_timestep,plot_in_sec)

    plt.tight_layout()
    if save_on:
        plt.savefig('%s/%s%s.png'%(save_dir,lab,fig_name),dpi=300)

def get_continuous_cmap(col_list,input_hex=False,float_list=None):
    ''' creates and returns a color map that can be used in heat map figures.
        If float_list is not provided, colour map graduates linearly between each color in col_list.
        If float_list is provided, each color in col_list is mapped to the respective location in float_list. 
        
        Parameters
        ----------
        col_list: list of color code strings
        float_list: list of floats between 0 and 1, same length as col_list. Must start with 0 and end with 1.
        
        Returns
        ----------
        colour map'''
    if input_hex:
        rgb_list = [rgb_to_dec(hex_to_rgb(i)) for i in col_list]
    else:
        rgb_list = col_list.copy()

    if float_list:
        pass
    else:
        float_list = list(np.linspace(0,1,len(rgb_list)))
        
    cdict = dict()
    for num, col in enumerate(['red', 'green', 'blue']):
        col_list = [[float_list[i], rgb_list[i][num], rgb_list[i][num]] for i in range(len(float_list))]
        cdict[col] = col_list
    cmp = mpl.colors.LinearSegmentedColormap('my_cmp', segmentdata=cdict, N=256)
    return cmp

def hex_to_rgb(value):
    '''
    Converts hex to rgb colours
    value: string of 6 characters representing a hex colour.
    Returns: list length 3 of RGB values'''
    value = value.strip("#") # removes hash symbol if present
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))

def rgb_to_dec(value):
    '''
    Converts rgb to decimal colours (i.e. divides each value by 256)
    value: list (length 3) of RGB values
    Returns: list (length 3) of decimal values'''
    return [v/256 for v in value]
