#!/usr/bin/env python3
'''
Functions related to plotting spatio-temporal evolution of variables as an image
By Jeena Yun
Last modification: 2025.10.08.
'''
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import cmcrameri.cm as cram
import myplots
mp = myplots.Figpref()

yr2sec = 365*24*60*60
mymint = (70/255,225/255,165/255)

# fields: Time [s] | state [s?] | cumslip0 [m] | traction0 [MPa] | slip-rate0 [m/s] | normal-stress [MPa]
# Index:     0     |      1     |       2      |        3        |         4        |          5 

def fout_image(outputs, dep, event_info, **kwargs):
    options = {
        'vmin': None,
        'vmax': None,
        'Vths': 1e-2,
        'zoom_frame': [],
        'horz_size': 20.6,
        'vert_size': 11,
        'sys_col' : mymint,
        'pr_col' : 'w',
        'no_classfication': False,
        'plot_in_timestep': True,
        'plot_in_sec': False,
        'colorbar_off': False,
        'scatter_off': False,
        'print_on': True,
        'save_on': True
    }
    options.update(kwargs)
    
    processed_outputs = class_figtype(outputs, event_info, options)[0]
    X, Y, var, options = get_var(processed_outputs, dep, options)
    ax = plot_image(X, Y, var, outputs, event_info, options)
    return ax

def get_var(outputs, dep, options):
    idx = {'state': 1, 'slip': 2, 'shearT': 3, 'delshearT': 3, 'sliprate': 4, 'normalT': 5, 'delnormalT': 5}
    ii = np.argsort(abs(dep))
    if options['image'] == 'dCFS':
        var = abs(outputs[ii,1:,3]) - 0.4 * abs(outputs[ii,1:,5])
        start_from_1 = True
    elif options['image'][-1] == 'T' or idx[options['image']] == 4:
        var = abs(outputs[ii,1:,idx[options['image']]])
        var = abs(var)
        start_from_1 = True
    else:
        var = outputs[ii,:,idx[options['image']]]
        start_from_1 = False

    if options['image'] == 'shearT': # ---- Total stress
        var = np.array([var[i,:] for i in range(var.shape[0])])
    elif 'd' in options['image']:    # ---- Stress change
        var = np.array([var[i,:]-var[i,0] for i in range(var.shape[0])])

    if options['plot_in_timestep']:
        print('Plot in time steps')
        xax = np.arange(var.shape[1])
    else:
        if options['plot_in_sec']:
            print('Plot in time [s]')
            xax = outputs[0,:,0]
        else:
            print('Plot in time [yr]')
            xax = np.array(outputs[0][:,0]/yr2sec)
        if start_from_1: 
            print('Sliprate or stresses - start from timestep 1 & take absolute')
            xax = xax[1:]

    X,Y = np.meshgrid(xax,np.sort(abs(dep)))
    options['start_from_1'] = start_from_1

    return X, Y, var, options

def gen_cmap(options):
    Vp = 1e-9
    V0 = 1e-6
    cb_dict = {'state': 'State Variable $\psi$', 
               'slip': 'Cumulative Slip [m]', 
               'shearT': 'Shear Traction [MPa]', 
               'delshearT': 'Shear Traction Change [MPa]', 
               'sliprate': 'Slip Rate [m/s]', 
               'normalT': 'Normal Stress [MPa]', 
               'delnormalT': 'Normal Stress Change [MPa]',
               'dCFS': '$\Delta$CFS [MPa]'} 
    if options['image'] == 'sliprate':
        cm = mpl.colormaps['RdYlBu_r']
        col_list = [cm(i)[0:3] for i in [0.15, 0.5, 0.8, 0.9]]
        col_list = [cm(0.15)[0:3], mpl.colormaps['jet'](0.67), cm(0.8)[0:3], cm(0.9)[0:3]]
        col_list.insert(0, (0, 0, 0))
        col_list.insert(3, mpl.colors.to_rgb('w'))
        col_list.append(mpl.colormaps['turbo'](0))
        float_list = [0, mpl.colors.LogNorm(options['vmin'], options['vmax'])(Vp), mpl.colors.LogNorm(options['vmin'], options['vmax'])(V0)]
        [float_list.append(k) for k in np.linspace(mpl.colors.LogNorm(options['vmin'], options['vmax'])(options['Vths']),1,4)]
        cmap_n = get_continuous_cmap(col_list, input_hex=False, float_list=float_list)
    elif 'd' in options['image']: # ---- plotting stress change
        cm = cram.vik
        col_list = [cm(i) for i in np.linspace(0,1,6)]
        col_list.insert(3,mpl.colors.to_rgb('w'))
        float_list = [mpl.colors.Normalize(options['vmin'], options['vmax'])(i) for i in np.linspace(options['vmin'], 0, 4)]
        [float_list.append(mpl.colors.Normalize(options['vmin'], options['vmax'])(k)) for k in np.linspace(0, options['vmax'], 4)[1:]]
        cmap_n = get_continuous_cmap(col_list, input_hex=False, float_list=float_list)
    elif options['image'][-1] == 'T' and 'del' not in options['image']: # ---- plotting total stress
        cmap_n = cram.davos
    elif options['image'] == 'state':
        cmap_n = 'magma'
    elif options['image'] == 'slip':
        cmap_n = 'viridis'
    cb_label = cb_dict[options['image']]
    return cmap_n,cb_label

def class_figtype(outputs, event_info, options):
    if len(event_info) == 0:
        options['scatter_off'] = True
        evdep, its_all, ite_all = [],[],[]
        if len(options['zoom_frame']) == 0: # Full image
            if options['print_on']: print('Full image')
            processed_outputs = outputs
            tsmin,tsmax,its,ite,buffer1,buffer2,iev1,iev2 = [],[],[],[],[],[],[],[]
            xl_opt = 1
            vlines_opt = 0
            xlab_opt = 1
            name_opt = 1
    else:
        tstart = event_info['tstart']
        tend = event_info['tend']
        evdep = event_info['evdep']
        its_all = np.array([np.argmin(abs(outputs[0][:,0]-t)) for t in tstart])
        ite_all = np.array([np.argmin(abs(outputs[0][:,0]-t)) for t in tend])
        if len(options['zoom_frame']) == 0: # Full image
            if options['print_on']: print('Full image')
            processed_outputs = outputs
            tsmin,tsmax,its,ite,buffer1,buffer2,iev1,iev2 = [],[],[],[],[],[],[],[]
            xl_opt = 1
            scatter_opt = 1
            vlines_opt = 0
            txt_opt = 0
            xlab_opt = 1
            name_opt = 1
    if options['scatter_off']: 
        scatter_opt, txt_opt = 0, 0
    fig_opts = [xl_opt,scatter_opt,vlines_opt,txt_opt,xlab_opt,name_opt]
    return processed_outputs,evdep,its_all,ite_all,tsmin,tsmax,its,ite,buffer1,buffer2,iev1,iev2,fig_opts

def decoration(outputs, event_info, Wf, acolor, options):
    # --- Figure type
    evdep,its_all,ite_all,tsmin,tsmax,its,ite,buffer1,buffer2,iev1,iev2,fig_opts = class_figtype(outputs, event_info, options)[1:]
    xl_opt,scatter_opt,vlines_opt,txt_opt,xlab_opt,name_opt = fig_opts

    if options['plot_in_sec']:
        if options['start_from_1']:
            time = np.array(outputs[0,1:,0])
        else:
            time = np.array(outputs[0,:,0])
    else:
        if options['start_from_1']:
            time = np.array(outputs[0,1:,0])/yr2sec
        else:
            time = np.array(outputs[0,:,0])/yr2sec

    if xl_opt == 1 or xl_opt == 2:
        xl = plt.gca().get_xlim()
    elif xl_opt == 3:
        width = np.round(time[ite]-time[its])
        xl = [time[its]-width/6, time[ite]+width/6]

    if scatter_opt == 2:
        lim_s, lim_e = tsmin, tsmax
    elif scatter_opt == 3:
        lim_s, lim_e = its, ite

    if len(its_all) > 0:
        if options['plot_in_timestep']:
            if scatter_opt == 1:
                xs = its_all
            else:
                xs = its_all - lim_s
        else:
            xs = time[its_all]

    len_event = 0
    if scatter_opt == 1:
        if len(evdep) > 0:
            plt.scatter(xs,evdep,s=700,marker='*',facecolor=options['sys_col'],edgecolor='k',lw=1.5,zorder=3,label='Earthquakes')
            len_event = len(evdep)
    elif scatter_opt > 0:
        ievents = np.where(np.logical_and(its_all>=lim_s,its_all<=lim_e))[0]
        plt.scatter(xs[ievents],evdep[ievents],s=700,marker='*',facecolor=options['sys_col'],edgecolor='k',lw=1.5,zorder=3,label='Earthquakes')
        len_event = len(ievents)

    if len_event > 0:
        plt.legend(fontsize=20,framealpha=1,loc='lower right')

    if vlines_opt == 1:
        if options['plot_in_timestep']:
            plt.vlines(x=[0,ite-its],ymin=0,ymax=Wf,linestyles='--',color=acolor,lw=1.5)
        else:
            plt.vlines(x=[time[its],time[ite]],ymin=0,ymax=Wf,linestyles='--',color=acolor,lw=1.5)
    elif vlines_opt == 2:
        if len(its_all) > 0:
            plt.vlines(x=xs[iev1:iev2+1],ymin=0,ymax=Wf,linestyles='--',color=acolor,lw=1.5)

    if len(its_all) > 0:
        if txt_opt == 1:
            for k in range(len(its_all)):
                if its_all[k]>=lim_s and its_all[k]<=lim_e:
                    plt.text(xs[k],evdep[k]-0.2,'%d'%(k),color=acolor, fontsize=20,ha='right',va='bottom')
        elif txt_opt == 2:
            plt.text(xs[iev1]+width/18,23,'Event %d'%(iev1), fontsize=30,fontweight='bold',color=acolor,ha='left',va='bottom') # coseismic
        elif txt_opt == 3:
            plt.text(xs[iev1]-width/150,12,'Event %d'%(iev1), fontsize=30,fontweight='bold',color=acolor,ha='right',va='bottom',rotation=90)
            if options['plot_in_timestep']:
                plt.text(ite-its+width/150,12,'Event %d'%(iev2), fontsize=30,fontweight='bold',color=acolor,ha='left',va='bottom',rotation=90)
            else:
                plt.text(time[ite]+width/150,12,'Event %d'%(iev2), fontsize=30,fontweight='bold',color=acolor,ha='left',va='bottom',rotation=90)

    if xl_opt == 1:
        if options['plot_in_timestep']:
            plt.xlim(0,xl[1])
        else:
            plt.xlim(np.min(time),np.max(time))
    elif xl_opt == 2:
        if options['plot_in_timestep']:
            plt.xlim(0,xl[1])
        else:
            plt.xlim(time[tsmin],xl[1])
    elif xl_opt == 3:
        plt.xlim(xl)
        
    plt.ylim(Wf, 0)
    plt.ylabel('Depth [km]', fontsize=30)

    if options['plot_in_timestep']:
        if xlab_opt == 1: 
            plt.xlabel('Timesteps', fontsize=30)
        elif xlab_opt == 2:
            if tsmin > 1e-3:
                xt = plt.gca().get_xticks()
                xt_lab = ['%d'%(ticks+tsmin) for ticks in xt]
                plt.gca().set_xticks(xt)
                plt.gca().set_xticklabels(xt_lab)
            plt.xlabel('Timesteps', fontsize=30)
    elif options['plot_in_sec']:        
        plt.xlabel('Time [s]', fontsize=30)
    else:        
        plt.xlabel('Time [yr]', fontsize=30)

    if name_opt == 1:
        fig_name = '_image'
    elif name_opt == 2:
        if np.mod(tsmin,tsmax-tsmin) == 0:
            fig_name = '_zoom_image_%d_%d'%(tsmin/(tsmax-tsmin),tsmax/(tsmax-tsmin))
        else:
            fig_name = '_zoom_image_%d_%d'%(tsmin/1000,tsmax/1000)
    elif name_opt == 3:
        fig_name = '_zoom_image_ev%d'%(iev1)
    elif name_opt == 4:
        fig_name = '_zoom_image_ev%dto%d'%(iev1,iev2)
    elif name_opt == 5:
        fig_name = '_zoom_image_interevent_%d_%d'%(iev1,iev2)

    return fig_name

def plot_image(X, Y, var, outputs, event_info, options):
    plt.rcParams['font.size'] = '27'
    _,ax = plt.subplots(figsize=(options['horz_size'], options['vert_size']))

    if 'd' in options['image']:
        maxdvar = max([abs(np.min(var)), abs(np.max(var))])
        options['vmin'] = -maxdvar
        options['vmax'] = maxdvar
    elif 'sliprate' in options['image']:
        options['vmin'] = 1e-12
        options['vmax'] = 1e1
    cmap_n,cb_label = gen_cmap(options)

    if options['image'] == 'sliprate':
        cb = plt.pcolormesh(X, Y, var, cmap=cmap_n, norm=mpl.colors.LogNorm(options['vmin'], options['vmax']))
        acolor = 'w'
    elif options['image'] == 'shearT': # ---- Total stress
        cb = plt.pcolormesh(X, Y, var, cmap=cmap_n, vmin= options['vmin'], vmax= options['vmax'])
        acolor = 'k'
    elif 'd' in options['image']:      # ---- Stress change
        cb = plt.pcolormesh(X, Y, var, cmap=cmap_n, vmin= options['vmin'], vmax= options['vmax'])
        acolor = 'k'
    else:
        cb = plt.pcolormesh(X, Y, var, cmap=cmap_n, vmin= options['vmin'], vmax= options['vmax'])
        acolor = 'w'
    if not options['colorbar_off']:
        plt.colorbar(cb, extend='both').set_label(cb_label, fontsize=30, rotation=270, labelpad=30) # vertical colorbar

    Wf = np.max(abs(Y))
    fig_name = decoration(outputs, event_info, Wf, acolor, options)

    plt.tight_layout()
    if  options['plot_in_timestep']:
        plt.savefig('%s/%s%s_timesteps.png'%(options['save_dir'], options['image'], fig_name),dpi=300)
    else:
        plt.savefig('%s/%s%s.png'%(options['save_dir'],options['image'], fig_name),dpi=300)
    return ax

def get_continuous_cmap(col_list, input_hex=False, float_list=None):
    import matplotlib.colors as mcolors
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
        float_list = list(np.linspace(0, 1, len(rgb_list)))
        
    cdict = dict()
    for num, col in enumerate(['red', 'green', 'blue']):
        col_list = [[float_list[i], rgb_list[i][num], rgb_list[i][num]] for i in range(len(float_list))]
        cdict[col] = col_list
    cmp = mcolors.LinearSegmentedColormap('my_cmp', segmentdata=cdict, N=256)
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