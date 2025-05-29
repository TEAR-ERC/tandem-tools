#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 15:39:10 2024

@author: bar
"""
import numpy as np
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import pandas as pd
#import gutenbergRichter
import cmap_tandem


                                    

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.ndimage import label, binary_dilation
#%%

class upperPlateFaultData:
    def __init__(self,data,functionForAValues=None):
        self.data=data
        ind=self.FindEQInd(np.abs(data['slip-rate0']))
        if functionForAValues is not None:
            self.a=functionForAValues(np.abs(data['z'].values))
        cfs=self.GetCFS()
        
        if len(ind)==0:
            self.avgcfs=self.ComputeAvgCFS(cfs[0:],data['Time'][0:])
        elif len(ind)==1:
            self.avgcfs=self.ComputeAvgCFS(cfs[0:ind[0]-100],data['Time'][0:ind[0]-100])
        elif len(ind)>1:
            self.avgcfs=self.ComputeAvgCFS(cfs[0:ind[0]-100],data['Time'][0:ind[0]-100])
        else:
            print(len(ind))
            raise ValueError("WTF")
        
        #self.avgcfs=365*24*365.25*self.avgcfs # cfs per year
        #self.avgcfs=self.ComputeAvgCFS(cfs[1000:5000])
        self.data=data.assign_coords(CFS=self.avgcfs)
        self.ind=ind

        
    def FindEQInd(self,slip_rate):
        ind=find_peaks(slip_rate,height=1e-5,distance=2000)[0]
        return ind
    
    def GetCFS(self):
        return self.data['normal-stress']*0.6+self.data['traction0']
    
   # def GetCFS(self):
    #    a=self.a
    #    fric=a*np.arcsinh((np.abs(self.data['slip-rate0'])/2*1e-6)*np.exp(self.data['state']/a))
    #    return fric*self.data['normal-stress']+self.data['traction0']
        
    
    def ComputeAvgCFS(self,cfs,time):
        return np.nanmean(np.gradient(cfs))
        #return (cfs[-1]-cfs[0])/(time[-1]-time[0])
        #return np.mean(np.cumsum(cfs))
    
    def PlotCFSAndSlip(self,ax=None):
        if ax is None:
            fig,ax=plt.subplots(2,1)
            
        cfs=self.GetCFS() 
        time=self.data['Time']
        slip_rate=np.log10(self.data['slip-rate0'])
        ax[0].plot(cfs)
        ax[1].plot(slip_rate)
        
        for ind_i in self.ind:
            #ax[0].scatter(time[ind_i],cfs[ind_i])
            #ax[1].scatter(time[ind_i],slip_rate[ind_i])
            x=np.arange(0,len(cfs) )
            ax[0].scatter(x[ind_i],cfs[ind_i])
            ax[1].scatter(x[ind_i],slip_rate[ind_i])
    def ReturnXrWithCFS(self):
        return self.data
    
#%%

class discerteFaults:
    def __init__(self,slipRate,slip,time,positionAlongFault,params=None):
        
        if slipRate.shape[0] != len(positionAlongFault):
            raise ValueError ("dimesntion of slip rate and fault distance do not align")
        if slipRate.shape[1] != len(time):
            raise ValueError ("dimesntion of slip rate and time do not align")
        if slipRate.shape != slip.shape:
            raise ValueError ("dimesntion of slip rate and slip do not align")
        if np.all(np.diff(positionAlongFault) >= 0) is False:
            raise ValueError ("position along fault must be sorted")

            
        self.slipRate=np.abs(np.asanyarray(slipRate))
        self.slip=np.abs(np.asanyarray(slip))
        self.time=np.asanyarray(time)
        self.positionAlongFault=np.asanyarray(positionAlongFault)
        self.params=params
        self.EQs=self.ComputeEQ()
        self.events=self.stickEventsTogther(self.EQs)
        self.ComputeCatalog()
       # try:
       #     self.GR=gutenbergRichter.gunterbergRichterFit(self.catalog['Mw'])  
       # except:
       #     print("cant fit GR")
                
        
    def Plot1d(self):
        x=np.arange(0, len(self.time))
        for j in range(len(self.slip[:,0])):
            EQs_j=self.EQs.loc[self.EQs['positionAlongFaultInd']==j]
            fig,ax=plt.subplots(1,2)
            ax[0].plot(x,self.slipRate[j,:])
            ax[1].plot(x,self.slip[j,:])
            ax[0].set_title(self.positionAlongFault[j])
            
            startInd=EQs_j['startInd'].values
            endInd=EQs_j['endInd'].values
            
            ax[0].scatter(x[startInd],self.slipRate[j,startInd],color='magenta')
            ax[0].scatter(x[startInd],self.slipRate[j,startInd],color='blue') 
            
            ax[1].scatter(x[endInd],self.slip[j,endInd],color='magenta')
            ax[1].scatter(x[endInd],self.slip[j,endInd],color='blue')
            
    def Plot2D(self):
        fig,ax=plt.subplots(1,1)
        cmap_coseismic=cmap_tandem.ReturnCmap(vmin=1e-19, vmax=1,Vths=1e-3)
        
        ax.imshow(np.log10(np.abs(self.slipRate)),extent=[0,len(self.time),np.max(self.positionAlongFault),0],aspect=1000,cmap=cmap_coseismic,vmin=-19,vmax=0)
        x=np.arange(0, len(self.time))
        for j in range(len(self.slip[:,0])):
            EQs_j=self.EQs.loc[self.EQs['positionAlongFaultInd']==j]
            
            startInd=EQs_j['startInd'].values
            endInd=EQs_j['endInd'].values
            
            posAlongFault=np.ones_like(startInd)*self.positionAlongFault[j]
            ax.scatter(startInd,posAlongFault,color='magenta')
            
    def Plot2DEvents(self):
        fig,ax=plt.subplots(1,1)
        cmap_coseismic=cmap_tandem.ReturnCmap(vmin=1e-19, vmax=1,Vths=1e-3)
        
        ax.imshow(np.log10(np.abs(self.slipRate)),extent=[0,len(self.time),np.max(self.positionAlongFault),0],aspect=1000,cmap=cmap_coseismic,vmin=-19,vmax=0)
        for event_i in self.events:
            ax.scatter(event_i['startInd'],event_i['positionAlongFault'])

        
    def ComputeEQ(self):
        EQs=[]
        for i in range(len(self.slipRate[:,0])):
            EQ_i=self.FindEQVector(self.slipRate[i,:],self.slip[i,:],i)
            if EQ_i is not None:
                
                EQs.append(EQ_i)
        try:        
            EQs=pd.concat(EQs,ignore_index=True)
        except ValueError:
            print("no EQs found")
            EQs = pd.DataFrame(columns=["startTime", "positionAlongFault", "positionAlongFaultInd"])

                
        return EQs
    
    def FindEQVector(self,slipRate,slip,positionAlongFaultInd):
    
        distance = self.params.get('distance', 10000)  # Default to 10000 if not provided
        height = self.params.get('height', 1e-4)       # Default to 1e-4 if not providedd
        EQ_duration = self.params.get('EQ_duration', 120)  # Default to 1 if not provided
            
        slipRatePeaks=find_peaks(np.gradient(np.gradient(slipRate,self.time),self.time),distance=distance,height=height)[0]
        
        if len(slipRatePeaks)==0:
            return None
            
        slipRaiseAlongFault=[] 
        for slipRatePeak_i in slipRatePeaks:
            onset=self.time[slipRatePeak_i]
            endOfEQ=np.argmin(np.abs((EQ_duration+onset-self.time[slipRatePeak_i:])))+slipRatePeak_i
            slip=self.slip[positionAlongFaultInd,endOfEQ]-self.slip[positionAlongFaultInd,slipRatePeak_i]
                #slipRaiseAlongFault.append(discerteFault(slipRate[slipRatePeak_i:endOfEQ],slip[slipRatePeak_i:endOfEQ],time=self.time[slipRatePeak_i:endOfEQ],positionAlongFault=positionAlongFault))
            slipRaiseAlongFault.append(pd.DataFrame(data={'startInd':slipRatePeak_i,'endInd':endOfEQ,'positionAlongFaultInd':positionAlongFaultInd,
                                                         'startTime':onset,'endTime':self.time[endOfEQ],
                                                         'positionAlongFault':self.positionAlongFault[positionAlongFaultInd],
                                                         'slip':slip},index=[0]))
            
        return pd.concat(slipRaiseAlongFault,ignore_index=True)
        
    def stickEventsTogther(self,EQs):
        events=[]
        
        EQs=EQs.sort_values(['startTime','positionAlongFault']).reset_index(drop=True)
        while len(EQs)>1:
            onset=EQs.loc[0,'startTime']
            tempTime=EQs.loc[:,'startTime']-onset-2000
            eventIndex = tempTime.loc[tempTime < 0].index
            newEvent=EQs.loc[eventIndex]
            newEvent=newEvent.sort_values('positionAlongFaultInd').reset_index(drop=True)
            events.append(newEvent)
            EQs.drop(eventIndex,inplace=True)
            EQs.reset_index(drop=True, inplace=True)
            
            
        return events
    
    
    def ComputeCatalog(self):
        alongDipEQExtent=[];L=[];Mw=[];startTime=[]
        for event in self.events:
            alongDipEQExtent_i=np.max(event.positionAlongFault)-np.min(event.positionAlongFault)
            if alongDipEQExtent_i==0:
                alongDipEQExtent_i=self.positionAlongFault[1]-self.positionAlongFault[0]
            integral = np.trapz(event.slip, event.positionAlongFault*1e3)
            L_i=(alongDipEQExtent_i/1.73)**1.412
            L.append(L_i)
            Mw.append(((np.log10(integral*L_i*1e3*30e9)-9.1)/1.5))
            startTime.append(np.min(event['startTime']))
            alongDipEQExtent.append(alongDipEQExtent_i)
            
        self.catalog=pd.DataFrame({'alongDipEQExtent':alongDipEQExtent,'onset':startTime,'alongStrikeExtent':L,'Mw':Mw})
            
            
               
               
               
               
        # def StickEQsTogther(self,EQs):
        #     EQs=EQs.flatten()
        #     indexForEqs=np.arange(0,len(EQs) )
            
        #     events=[]
            
        #     for i in indexForEqs:
        #         for j in indexForEqs:
        #             if i==j or EQs[i].faultNode == EQs[j].faultNode:
        #                 break
                    
        #             EQ1_onset=self.time[EQs[i].start]
        #             EQ2_onset=self.time[EQs[j].start]
                    
        #             if np.abs(EQ1_onset-EQ2_onset) <2000:
        #                 for event_i in events:
        #                     for EQ_l_in_event in event_i:
        #                         if (np.abs(EQ1_onset-EQ_l_in_event) <2000) or (np.abs(EQ2_onset-EQ_l_in_event) <2000):
        #                             event_i.append(EQs[i])
        #                             event_i.append(EQs[j])
                                    
        #                             foundMatchToEvent=True
        #                             break;
                                    
        #                     if foundMatchToEvent:
        #                         break
                            
        #                 events.append(EQs[i],EQs[j])
                        
        #                 EQs.drop(i,j)
                                    
                                    
                                    
                                    


class ConnectedComponentEQ:
    """
    Earthquake detector based on 2‑D (depth × time) connected‑component labelling.

    Parameters
    ----------
    slip_rate : 2‑D array  [n_dip , n_t]
    slip      : 2‑D array  same shape (cumulative slip)
    time      : 1‑D array  length n_t   (s, yr — arbitrary)
    dip_pos   : 1‑D array  distance/depth along dip  (km or m, monotonic)
    opts      : dict with optional keys
        rate_threshold : float   (default 1e‑3, units of slip_rate)
        connectivity   : {1,2}   1=4‑neighbour, 2=8‑neighbour   (default 2)
        mu             : float   shear modulus (Pa, default 30e9)
        min_pixels     : int     discard blobs smaller than this (default 3)
    """

    def __init__(self, slip_rate, slip, time, dip_pos, *, opts=None):
        # --------------- store & sanity‑check ----------------------
        if slip_rate.shape != slip.shape:
            raise ValueError("slip_rate and slip must have identical shape")
        if len(dip_pos) != slip_rate.shape[0]:
            raise ValueError("dip_pos length ≠ first slip_rate dim")
        if len(time) != slip_rate.shape[1]:
            raise ValueError("time length ≠ second slip_rate dim")
        if np.any(np.diff(dip_pos) < 0):
            raise ValueError("dip_pos must be monotonically increasing")

        self.V   = np.abs(np.asarray(slip_rate))
        self.S   = np.abs(np.asarray(slip))
        self.t   = np.asarray(time)
        self.dip = np.asarray(dip_pos)

        o = opts or {}
        self.Vthr  = o.get("rate_threshold", 1e-3)
        self.conn  = o.get("connectivity", 2)
        self.mu    = o.get("mu", 30e9)
        self.minpx = o.get("min_pixels", 3)

        # --------------- main pipeline ----------------------------
        self.mask, self.labels, self.nlbl = self._build_labels()
        self.catalog = self._build_catalog()

    # ==============================================================
    # internal helpers
    # ==============================================================

    def _build_labels(self):
        """Binary mask + connected components."""
        mask = np.abs(self.V) > self.Vthr
        struct = np.ones((3,3), bool) if self.conn == 2 else \
                 np.array([[0,1,0],[1,1,1],[0,1,0]], bool)
        labels, n = label(mask, structure=struct)

        # size filter
        if self.minpx > 1:
            sizes = np.bincount(labels.ravel())
            bad   = np.where(sizes < self.minpx)[0]
            for b in bad:
                labels[labels == b] = 0
            # relabel to 1 … n_clean
            labels, n = label(labels > 0, structure=struct)

        return mask, labels, n

    def _slip_for_blob(self, blob_mask):
        """
        Average coseismic slip for one event.
        For each dip index, use slip difference between first & last
        time sample that belong to the blob.
        """
        dip_idx, time_idx = np.where(blob_mask)
        slip_depth = {}
        for d, t in zip(dip_idx, time_idx):
            slip_depth.setdefault(d, []).append(t)
        slip_inc = []
        for d, times in slip_depth.items():
            t0, t1 = min(times), max(times)
            slip_inc.append(self.S[d, t1] - self.S[d, t0])
        return np.mean(slip_inc)

    def _build_catalog(self):
        """Return a Pandas DF with Mw etc."""
        rows = []
        for k in range(1, self.nlbl+1):
            m = self.labels == k
            dip_idx, time_idx = np.where(m)
            dmin, dmax = self.dip[dip_idx].min(), self.dip[dip_idx].max()
            along_dip  = dmax - dmin                       # same units as dip_pos
            if along_dip == 0:                             # single row ⇒ use grid spacing
                along_dip = np.diff(self.dip).mean()
            L = (along_dip / 1.73) ** 1.412                # along‑strike length
            slip_bar = self._slip_for_blob(m)
            area = along_dip * L * 1e6                     # km² → m²  (1e3×1e3)
            M0 = self.mu * slip_bar * area                 # N·m
            Mw = (np.log10(M0) - 9.1) / 1.5
            rows.append(dict(label=k,
                             t_start=self.t[time_idx.min()],
                             t_end  =self.t[time_idx.max()],
                             depth_top=dmin, depth_bot=dmax,
                             along_dip=along_dip, along_strike=L,
                             slip_mean=slip_bar,
                             Mw=Mw))
        return pd.DataFrame(rows).sort_values("t_start").reset_index(drop=True)

    # ==============================================================
    # public helpers
    # ==============================================================

    def plot_events(self, cmap_base='viridis', alpha_mask=0.4, figsize=(10,6)):
        """
        Show log10|slip_rate| and outline each detected EQ in a unique colour.
        """
        fig, ax = plt.subplots(figsize=figsize)

        im = ax.imshow(np.log10(np.abs(self.V)),
                       extent=[self.t[0], self.t[-1], self.dip[-1], self.dip[0]],
                       aspect='auto', cmap=cmap_base)
        ax.set_xlabel("Time")
        ax.set_ylabel("Distance along dip (km)")
        cb = fig.colorbar(im, ax=ax, label='log₁₀|slip‑rate|')

        # coloured perimeters
        colours = cm.get_cmap('tab10', self.nlbl)
        for k in range(1, self.nlbl+1):
            m = self.labels == k
            edge = binary_dilation(m) ^ m           # simple perimeter
            dip_idx, time_idx = np.where(edge)
            ax.scatter(self.t[time_idx], self.dip[dip_idx],
                       s=4, c=[colours(k-1)], marker='s', linewidths=0)

        ax.set_title("Detected coseismic ruptures")
        fig.tight_layout()
        return fig, ax
                              
                        
            
            
            
        
        

        
            
            
            