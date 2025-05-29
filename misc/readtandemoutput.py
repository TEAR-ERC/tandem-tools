#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 15:35:32 2024

@author: bar
"""

import pandas as pd
import xarray as xr
import numpy as np
import glob
import matplotlib.pyplot as plt
import os
import netCDF4
import discerteFaultsTandem as dsft
#import gutenbergRichter
import cmap_tandem
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.gridspec import GridSpec
#%%

    
def load_csv_extract_x_position(file_path):
    """
    Reads the first line of a CSV file to extract the position of 'x',
    then loads the rest of the CSV file into a pandas DataFrame, skipping the first line.

    Parameters:
    - file_path: str, path to the CSV file.

    Returns:
    - x_position: list of floats, the extracted position of 'x'.
    - df: pandas DataFrame, the loaded CSV data excluding the first line.
    """

    # Open the file and read the first line
    with open(file_path, 'r') as file:
        first_line = file.readline().strip()

    # Extract the position of 'x' using a regular expression
    #match = re.search(r'x = \[([-\d.e]+), ([-\d.e]+)\]', first_line)
    #if match:
    #    x_position = [float(match.group(1)), float(match.group(2))]
    #else:
    #    x_position = None
    #    print("x position not found in the first line.")
    position=first_line.split(']')[0].split('[')[-1].split(',')


    x=float(position[0])
    y=float(position[1])

    if len(position)==2:
        xyz=np.array([x,y])
    elif len(position)==3:
        z=float(position[2])
        xyz=np.array([x,y,z])
    else:
        raise ("problem  extracting position data from csv file header")

    # Load the rest of the CSV file into a pandas DataFrame, skipping the first line
    df = pd.read_csv(file_path, skiprows=1)

    return xyz, df
#%%

def pandas_to_xarray(df, x_coord,z_coord, y_coord=None):
    """
    Reads a CSV file with a 'time' column and multiple data columns, and converts it into an xarray Dataset.

    Parameters:
    - filepath: str, path to the CSV file.
    - x_coord: float, x coordinate of the location.
    - y_coord: float, y coordinate of the location.

    Returns:
    - ds: xarray Dataset containing the data from the CSV file.
    """


    # Ensure 'time' column is a datetime type
    #df['time'] = pd.to_datetime(df['time'])

    # Initialize an empty list to store the DataArrays
    data_arrays = []

    # Iterate over the columns in the DataFrame, skipping the 'time' column
    for column in df.columns[1:]:  # Assuming the first column is 'time'
        # Create a DataArray for the current column
        if y_coord is None:
            da = xr.DataArray(df[column].values, dims=['Time'],
                        coords={'Time': df['Time'],
                                'x': x_coord,
                                'z': z_coord},
                        name=column)
        else:
            da = xr.DataArray(df[column].values, dims=['Time'],
                        coords={'Time': df['Time'],
                                'x': x_coord,
                                'z': z_coord,'y':y_coord},
                        name=column)
        data_arrays.append(da)

    # Combine the DataArrays into a single Dataset
    ds = xr.merge(data_arrays)
    return ds
#%%
def ReturnDataSets(path,pattern="fltst_*"):
    files=sorted(glob.glob(path+pattern))
    ds=[]
    for file in files:
        position,df=load_csv_extract_x_position(file)
        if len(position)==2:
            ds.append(pandas_to_xarray(df,position[0],position[1]))
        else:
            ds.append(pandas_to_xarray(df,position[0],position[1],position[2]))



    return ds

#%%
def PlotSlipFrequnelty(slip,time,frequency,z,plotAttributes=None,ax=None,coseismic=False):
    if ax is None:
        fig,ax=plt.subplots()

    i=0
    while i<(len(time)-1):

        if coseismic is True:
            if time[i+1]-time[i]>1:
                i+=1
                continue
        j=i+1
        while (time[j]-time[i]<frequency) & (j<len(time)-1):
            j+=1

        ax.plot(slip[:,j],z,**plotAttributes)
        i=j+1


#%%
def ReturnManyFaults(path,pattern='fltst_*',sperator="-",dim='z'):
    files=sorted(glob.glob(path+"*"))
    file_names = [os.path.basename(file) for file in files]
    fileheader=[]
    faults_data=[]

    if sperator is None:  ##just one fault
        data=ReturnDataSets(path=path,pattern=pattern)
        faults_data.append(faultObject(xr.concat(data, dim=(dim))))
        return faults_data


    for file_i in file_names:
        fileheader.append(file_i.split("-")[0])

    faults=np.unique(fileheader)


    for fault_i in faults:
        data=ReturnDataSets(path=path,pattern=fault_i+"*")
        faults_data.append(faultObject(xr.concat(data, dim=(dim))))

    return faults_data
#%%
def ReturnDomainProbeObject(path,pattern='dotst_domain*'):
    csv_datasets=ReturnDataSets(path,pattern)
    domainProbeDataset=xr.concat(csv_datasets, dim='x')
    return domainProbe(domainProbeDataset)
#%%
def FaultWork(slip, traction, distance):
    """Calculates the integral of traction * slip over distance using the trapezoidal rule."""
    return np.trapz(traction * slip, x=distance)

#%%
class outputBasedXarray:
    def __init__(self, dataset: xr.Dataset):
        """Initialize the Fault class with a given xarray.Dataset."""
        if not isinstance(dataset, xr.Dataset):
            raise TypeError(f"Expected dataset to be of type xr.Dataset, got {type(dataset)} instead.")
        self.dataset = dataset

    def SaveFaultData(self,filename,complevel=9,dtype="float32"):
        compression_settings = {var: {"zlib": True, "complevel": complevel,"shuffle": True, "dtype": dtype} for var in self.dataset.data_vars}
        if os.path.exists(filename):
            os.remove(filename)

        self.dataset.to_netcdf(filename, encoding=compression_settings, engine="netcdf4", mode="w")
#%%
class domainProbe(outputBasedXarray):
    def PlotDataset(self,ax=None,quantity='u1',plottingProperties=None,ratio=1):
        if ax is None:
            fig, ax = plt.subplots()

        if plottingProperties is None:
            zmin=np.min(self.dataset['x'].values);zmax=np.max(self.dataset['x'].values);y_range = np.abs(zmax - zmin)
            maxTime=len(self.dataset['Time'].values)
            extent=[0,maxTime,zmin,zmax];aspect = maxTime / y_range
            plottingProperties={'aspect':ratio*aspect,'extent':extent}


        dataToPlot=self.dataset[quantity]


        cb=ax.imshow(dataToPlot.values,**plottingProperties)
        return cb

    def Plot2DData(self,ax=None,quantity='u1',plottingProperties=None,time_in_yrs=None):
        if ax is None:
            fig, ax = plt.subplots()
        if time_in_yrs is None:
            nearest_index=-1
        else:
            nearest_index = abs(self.dataset['Time']//3600/24/365.25 - time_in_yrs).argmin().item()

        if plottingProperties is None:
            plottingProperties={}

        valueToPlot=self.dataset[quantity].values[:,nearest_index]



        x=self.dataset['x']
        ax.plot(x,valueToPlot,**plottingProperties)



#%%
class faultObject(outputBasedXarray):
    def __init__(self, dataset: xr.Dataset, a_func=None, b_func=None):
        """Initialize the Fault class with a given xarray.Dataset."""
        if not isinstance(dataset, xr.Dataset):
            raise TypeError(f"Expected dataset to be of type xr.Dataset, got {type(dataset)} instead.")
        super().__init__(dataset)
        self.RemoveDuplicates()
        self.RemoveNonMonotonicIndices()
        self._validate_dataset()
        self.distance_along_fault = self._compute_distance()

        #self.add_shear_integral()
        # Use the provided functions or fallback to defaults.
        self.a_func = a_func if a_func is not None else self.default_a_func
        self.b_func = b_func if b_func is not None else self.default_b_func
        self.dataset['CMF']=self.compute_failure_residual(0.005699833)
        self.dataset['staticCMF']=self.ComputeStaticCMF(0.6)

    def default_a_func(self, x, y):
        """
        Default function for friction parameter a(x,y).

        Uses the absolute value of y as depth (d) and defines:
            a0   = 0.010
            amax = 0.025
            H0   = 4.0
            H    = 8.0
            h    = 6.0

        The function returns:
            - a = a0 + (amax - a0)*(H0 - d)/H0    for d < H0,
            - a = a0                              for H0 <= d < H,
            - a = a0 + (amax - a0)*(d - H)/h         for H <= d < H+h,
            - a = amax                            for d >= H+h.
        """
        a0   = 0.010
        amax = 0.025
        H0   = 4.0
        H    = 8.0
        h    = 6.0

        d = abs(y)
        if d < H0:
            return a0 + (amax - a0) * (H0 - d) / H0
        elif d < H:
            return a0
        elif d < H + h:
            return a0 + (amax - a0) * (d - H) / h
        else:
            return amax

    def default_b_func(self, x, y):
        """
        Default function for friction parameter b(x,y).
        This function simply returns a constant value of 0.015.
        """
        return 0.015


    def compute_failure_residual(self,  D_RS, V0=1e-6 ,f0=0.6):
        """
        Compute the residual defined as:

            residual = τ - σₙ * F(V,θ)

        where:

            F(V, θ) = a(x,y) * asinh{ [|V|/(2V0)] * exp((f0 + b(x,y)*ln(V0*θ/D_RS)) / a(x,y) ) }

        Each row in the provided 2D arrays corresponds to a unique (x,y) position,
        and columns correspond to time.

        Parameters:
            V0 (float): Reference slip rate.
            D_RS (float): Effective critical slip distance.
            f0 (float): Reference friction value (default is 0.6).

        Returns:
            xr.DataArray: An array (with the same coordinates and dimensions as the input "state" variable)
                          containing the computed residual.
        """
        # Extract required 2D arrays from the dataset.
        # These arrays are assumed to have shape (n_positions, n_time).
        state         = self.dataset["state"].values         # θ, shape (n_positions, n_time)
        slip_rate     = self.dataset["slip-rate0"].values     # V, shape (n_positions, n_time)
        traction      = self.dataset["traction0"].values       # τ, shape (n_positions, n_time)
        normal_stress = self.dataset["normal-stress"].values   # σₙ, shape (n_positions, n_time)

        # Get spatial coordinates (each row corresponds to a unique (x,y) position).
        x_coords = self.dataset.coords["x"].values  # shape (n_positions,)
        y_coords = self.dataset.coords["z"].values  # shape (n_positions,)

        # Initialize the residual array.
        residual = np.empty_like(state)

        # Loop over each fault element (each row).
        for i in range(state.shape[0]):
            # Compute spatial friction parameters for this element.
            a_val = self.a_func(x_coords[i], y_coords[i])
            b_val = self.b_func(x_coords[i], y_coords[i])

            # Get the time series for the current fault element.
            V_i     = slip_rate[i, :]   # slip rate time series
            theta_i = state[i, :]       # state variable time series

            # Compute F(V,θ) for each time step.
            F_val = a_val * np.arcsinh(
                (np.abs(V_i) / (2 * V0)) *
                np.exp((f0 + b_val * np.log(V0 * theta_i / D_RS)) / a_val)
            )

            # Compute friction threshold: σₙ * F(V,θ)
            friction_threshold = normal_stress[i, :] * F_val

            # Because traction is negative for compression, we use -traction to represent its magnitude.
            # Compute the residual: (-τ) - (σₙ * F(V,θ))
            residual[i, :] = (-traction[i, :]) - friction_threshold

        # Wrap the result in an xarray.DataArray with the same coordinates and dims as "state"
        residual_da = xr.DataArray(residual, coords=self.dataset["state"].coords, dims=self.dataset["state"].dims)
        return residual_da


    def ComputeStaticCMF(self,f0=0.6):
        return  -self.dataset["traction0"] - self.dataset["normal-stress"]*f0


    def InterpolateToNewFaultObject(self,N=1000):
        z_new = np.linspace(self.dataset.z.min().item(), self.dataset.z.max().item(), N)
        dataset_interpolated = self.dataset.interp(z=z_new)
        return faultObject(dataset_interpolated)


    def _validate_dataset(self):
        """Check if required coordinates and data variables are present."""
        missing = {'x', 'z', 'slip0', 'traction0'} - set(self.dataset.coords) - set(self.dataset.data_vars)
        if missing:
            raise ValueError(f"Dataset is missing: {missing}")

        is_increasing = (self.dataset['Time'].diff(dim='Time') > 0).all()
        #if ~is_increasing:
        #    raise ValueError("Time coordinate is not monotonic.")

    def _compute_distance(self):
        """Calculate cumulative distance along the fault based on x and z coordinates."""
        x, z = self.dataset['x'].values, self.dataset['z'].values
        return np.concatenate(([0], np.cumsum(np.sqrt(np.diff(x)**2 + np.diff(z)**2))))

    def calculate_shear_integral(self):
        """Compute the shear integral (|traction0| * slip0 * dx) for each time step."""
        slip, traction = self.dataset['slip0'], np.abs(self.dataset['traction0'] * 1e6)
        return xr.apply_ufunc(
            lambda s, t: FaultWork(s, t, self.distance_along_fault*1e3),
            slip, traction, input_core_dims=[['z'], ['z']], vectorize=True
        ).rename('shear_integral')

    def ComputeAverageWork(self):
        """ Compute the average rate of work (shear_integral per unit time per unit distance) by integrating over time and diving by total time and length.
        Returns: The average work per unit time & distance.    """

        total_work = self.dataset['shear_integral'].integrate('Time')
        time_span = self.dataset['Time'].max() - self.dataset['Time'].min()
        return ((total_work/time_span).item())/(np.max(self.distance_along_fault)-np.min(self.distance_along_fault))

    def ComputeCoulombStressChange(self,fric):
        self.dataset['CFS']=self.dataset['traction0']+fric*self.dataset['normal-stress']

    def add_shear_integral(self):
        """Add shear_integral to the dataset."""
        self.dataset['shear_integral'] = self.calculate_shear_integral()
        return self.dataset



    def plot_shear_integral(self, ax=None):
        """Plot shear integral over time using the provided axis, or create one if None."""
        if ax is None:
            fig, ax = plt.subplots()

        self.dataset['shear_integral'].plot.line(ax=ax, x='Time', marker='.', linestyle='-')
        ax.set_xlabel('Time')
        ax.set_ylabel('Shear Integral (Nm)')
        ax.set_title('Shear Integral Over Time')
        ax.grid(True)

    def SelectPartOfData(self,min_time_in_yrs,max_time_in_yrs):
        if min_time_in_yrs >= max_time_in_yrs:
            raise ValueError("min_time_in_yrs must be less than max_time_in_yrs.")
        seconds_per_year = 365.25 * 24 * 3600
        min_time = min_time_in_yrs * seconds_per_year
        max_time = max_time_in_yrs * seconds_per_year

        mask_min = self.dataset['Time'] >= min_time
        mask_max = self.dataset['Time'] <= max_time

        combined_mask = mask_min & mask_max
        ds_subset = self.dataset.where(combined_mask, drop=True)
        #try :
        #    ds_subset = self.dataset.sel(Time=slice(min_time, max_time))
        #except KeyError:
        #    min_time_nearest = self.dataset['Time'].sel(Time=min_time, method="nearest").item()
        #    max_time_nearest = self.dataset['Time'].sel(Time=max_time, method="nearest").item()
        #    ds_subset = self.dataset.sel(Time=slice(min_time_nearest, max_time_nearest))

        return (faultObject(ds_subset))

    def ComputeCatalog(self,paramsToComputeCatalog=None):
        if paramsToComputeCatalog is None:
            paramsToComputeCatalog={"height":1e-3,'distance':1500}

        self.catalogObject=dsft.discerteFaults(self.dataset['slip-rate0'].values,self.dataset['slip0'].values,self.dataset['Time'].values,self.distance_along_fault,params=paramsToComputeCatalog)
        self.catalog=self.catalogObject.catalog


    def PlotDataset(self,ax=None,quantity='slip-rate0',plottingProperties=None,ratio=1):
        if ax is None:
            fig, ax = plt.subplots()
        creatPlottingProp=True
        if plottingProperties is not None:
            creatPlottingProp=False

        if creatPlottingProp:
            zmin=np.min(self.dataset['z'].values);zmax=np.max(self.dataset['z'].values);y_range = np.abs(zmax - zmin)
            maxTime=len(self.dataset['Time'].values)
            extent=[0,maxTime,zmin,zmax];aspect = maxTime / y_range
            plottingProperties={'aspect':ratio*aspect,'extent':extent}


        dataToPlot=self.dataset[quantity]

        if quantity=='slip-rate0' and creatPlottingProp:
            dataToPlot=np.log10(np.abs(dataToPlot))
            cmap_coseismic=cmap_tandem.ReturnCmap(vmin=1e-18, vmax=1,Vths=1e-3)
            plottingProperties['cmap']=cmap_coseismic
            plottingProperties['vmin']=-18
            plottingProperties['vmax']=0

        cb=ax.imshow(dataToPlot.values,**plottingProperties)
        return cb

    def GetRcov(self):
        R=np.diff(self.catalog.onset)/365.25/24/3600
        Rmean=np.mean(R)
        Rcov=np.std(R)

        return Rcov/Rmean


    def Plot2DData(self,ax=None,quantity='slip-rate0',plottingProperties=None,z=None,time_in_yrs=True):
        if ax is None:
            fig, ax = plt.subplots()
        if z is None:
            nearest_index=0
        else:
            nearest_index = abs(self.dataset.z - z).argmin().item()

        if plottingProperties is None:
            plottingProperties={}

        valueToPlot=self.dataset[quantity].values[nearest_index,:]

        if quantity=='slip-rate0':
            valueToPlot=np.log10(np.abs(valueToPlot))

        if time_in_yrs:
            x=self.dataset['Time']/3600/24/365.25
        else:
            x=np.arange(len(valueToPlot) )
        ax.plot(x,valueToPlot,**plottingProperties)


    def PlotEventsCatalogFault(self,ax=None,plottingProperties=None,ymin=None,ymax=None,timeSteps=False):
        if ymin is None and ymax is None and ax is not None:
            ymin, ymax = ax.get_ylim()
        if ax is None:
            fig, ax = plt.subplots()

        if plottingProperties is None:
            plottingProperties={}
        if timeSteps:
            timeIndex=self.Get_closest_time_index(self.catalog['onset'].values/365.25/3600/24)
            ax.vlines(timeIndex, ymin=ymin, ymax=ymax,**plottingProperties )
        else:
            ax.vlines(self.catalog['onset']/365.25/3600/24, ymin=ymin, ymax=ymax,**plottingProperties )


    def PlotMaxVel(self,ax=None,plottingProperties=None,time_in_yrs=True):
        if ax is None:
            fig, ax = plt.subplots()
        if plottingProperties is None:
            plottingProperties={}

        maxVel=np.log10(np.max(self.dataset['slip-rate0'].values,axis=0))
        if time_in_yrs:
            x=self.dataset['Time']/3600/24/365.25
        else:
            x=np.arange(len(maxVel))
        ax.plot(x,maxVel,**plottingProperties)




    def RemoveDuplicates(self):
        """
        Remove duplicates from self.dataset['Time'] and print the number of duplicates removed.
        """
        time_values = self.dataset['Time'].to_index()
        duplicates_mask = time_values.duplicated(keep='first')
        num_duplicates = duplicates_mask.sum()

        if num_duplicates > 0:
            #self.dataset = self.dataset.sel(Time=~duplicates_mask)
            print(f"Found {num_duplicates} duplicated Time values removed.")



    def RemoveNonMonotonicIndices(self):
        """
        Remove non-monotonic Time indices from the dataset and print the number removed.
        Assumes Time should be strictly increasing.
        """
        time = self.dataset['Time'].values
        # Create a mask where each Time is greater than the previous one
        mask = time == np.maximum.accumulate(time)
        num_non_monotonic = (~mask).sum()
        if num_non_monotonic > 0:
            #self.dataset = self.dataset.isel(Time=mask)
            print(f"Found {num_non_monotonic} non-monotonic Time values.")

    def ReturnQuantity(self,quantity='slip-rate0',z=0):
        if z is None:
            nearest_index=0
        else:
            nearest_index = abs(self.dataset.z - z).argmin().item()

        return self.dataset[quantity].values[nearest_index,:]



    def Get_closest_time_index(self, target_times):
        """
        Find the indices of the closest Time values to the given target times (in years).
        Uses absolute differences and argmin for each target time.
        """
        oneIndex=False
        if not isinstance(target_times, (list, np.ndarray)):
            target_times = [target_times]
            oneIndex=True

        # Convert target_times from years to seconds
        target_times_in_seconds = [t * 365.25 * 24 * 3600 for t in target_times]

        closest_indices = []


        for t in target_times_in_seconds:

            # Calculate absolute differences between each Time value and the target time
            time_diff = abs(self.dataset['Time'] - t)

            # Find the index of the smallest difference
            closest_index = time_diff.argmin().item()

            # Append the closest index to the results list
            closest_indices.append(closest_index)

        if oneIndex:
            closest_indices=closest_indices[0]

        return closest_indices











#%%
class domain:
    def __init__(self, faults,colors=None):
        self.faults=faults
        self.colors=colors
        #self.sort_faults_by_xshallow()
        #self.ComputeColors()


    def sort_faults_by_xshallow(self):
        def get_xshallow(fault):
            ID = np.argmax(fault.dataset.z.values)
            return fault.dataset.x.values[ID]

        self.faults.sort(key=get_xshallow)


    def IntepolateToHigherResoultion(self,incraseByHowMuch=5):
        intepolatedFaults=[]
        for fault_i in self.faults:
            N=len(fault_i.dataset.z)
            iterpolatedFault=fault_i.InterpolateToNewFaultObject(N*incraseByHowMuch)
            intepolatedFaults.append(iterpolatedFault)

        return intepolatedFaults

    def ComputeColors(self,interfaceColor='r'):
        avgwork=[]
        for i,fault_i in enumerate(self.faults):
            avgwork.append(fault_i.ComputeAverageWork())

        if len(self.faults)>2:
            norm = Normalize(vmin=np.min(np.log10(avgwork[1:])), vmax=np.max(np.log10(avgwork[1:])))
            cmap = plt.cm.plasma
            self.colors = cmap(norm(np.log10(avgwork)))
            self.colors[0,:]=np.array([1, 0, 0, 1])
        elif len(self.faults)==2:
            norm = Normalize(vmin=np.min(avgwork), vmax=np.max(avgwork))
            cmap = plt.cm.plasma
            self.colors = cmap(norm(avgwork))
            self.colors[0,:]=np.array([1, 0, 0, 1])
        else:
            self.colors=np.array([1, 0, 0, 1]).reshape(1, 4)

    def ComputeCatalog(self,colors,paramsToComputeCatalog=None):
        Mw=[]
        colorsPerEQs=[]
        for i,fault_i in enumerate(self.faults):
            fault_i.ComputeCatalog(paramsToComputeCatalog)
            fault_i.catalogObject.Plot2DEvents()
            filtered_values = fault_i.catalog['Mw'].copy()[np.isfinite(fault_i.catalog['Mw'])]
            Mw.append(filtered_values)
            if isinstance(colors, list):
                repeated_color = [colors[i]] * len(fault_i.catalog['Mw'])
            else:
                raise TypeError("colors must be a list.")
            colorsPerEQs.append(repeated_color)
            #colorsPerEQs.append(np.repeat(colors[i, :][np.newaxis, :], len(fault_i.catalog['Mw']), axis=0))

        self.Mw = pd.concat(Mw, ignore_index=True)
        #self.GR=gutenbergRichter.gunterbergRichterFit(self.Mw)

        self.colorsPerEQs= [color for sublist in colorsPerEQs for color in sublist]

    def SelectPartOfData(self,min_time_in_yrs,max_time_in_yrs):
        newFaults=[]
        for fault_i in self.faults:
            newFaults.append(fault_i.SelectPartOfData(min_time_in_yrs,max_time_in_yrs))

        return domain(newFaults,colors=self.colors)

    def PlotQuantity(self,ratio=0.5,quantity='slip-rate0'):

        fig,ax=plt.subplots(len(self.faults))


        cmap_coseismic=cmap_tandem.ReturnCmap(vmin=1e-19, vmax=1,Vths=1e-3)
        maxTime=len(self.faults[0].dataset['Time'].values)

        if len(self.faults)==1:
            ax=[ax]

        for i,fault_i in enumerate(self.faults):

            zmin=np.min(fault_i.dataset['z'].values);zmax=np.max(fault_i.dataset['z'].values);y_range = np.abs(zmax - zmin)
            extent=[0,maxTime,zmin,zmax];aspect = maxTime / y_range
            #plottingProperties={'vmin':-19,'vmax':0,'aspect':ratio*aspect,'extent':extent,'cmap':cmap_coseismic}
            #cb=fault_i.PlotDataset(ax=ax[i],quantity='slip-rate0',plottingProperties=plottingProperties)
            fault_i.PlotDataset(quantity=quantity,ax=ax[i])
        return fig


    def PlotQuantityNew(self, maxDepth=20, quantity='normal-stress', extraPlotPror=None):
    # 1) compute global vmin/vmax
        all_vals = np.hstack([f.dataset[quantity].values.ravel() for f in self.faults])
        global_vmin, global_vmax = all_vals.min(), all_vals.max()

        # detect if a Norm was provided
        use_norm = bool(extraPlotPror and 'norm' in extraPlotPror)

        # 2) override only if not using Norm
        vmin = extraPlotPror.get('vmin', global_vmin) if extraPlotPror else global_vmin
        vmax = extraPlotPror.get('vmax', global_vmax) if extraPlotPror else global_vmax

        # compute heights & time
        maxTime = len(self.faults[0].dataset['Time'].values)
        heights = []
        for f in self.faults:
            z = f.dataset['z'].values
            zmin, zmax = z.min(), z.max()
            if abs(zmin) > abs(maxDepth):
                zmin = -abs(maxDepth)
            heights.append(abs(zmax - zmin))

        fig = plt.figure(figsize=(12, sum(heights)))
        gs  = GridSpec(len(self.faults), 1, height_ratios=heights)
        axes, ims = [], []

        for i, fault in enumerate(self.faults):
            z = fault.dataset['z'].values
            zmin, zmax = z.min(), z.max()
            if abs(zmin) > abs(maxDepth):
                zmin = -abs(maxDepth)
            extent = [0, maxTime, zmin, zmax]

            ax = fig.add_subplot(gs[i, 0])
            axes.append(ax)

            # base kwargs
            p = {'aspect':'auto', 'extent':extent}

            # if no Norm, inject vmin/vmax
            if not use_norm:
                p.update(vmin=vmin, vmax=vmax)

            # merge any extras (this will include norm if present)
            if extraPlotPror:
                p.update(extraPlotPror)

            im = ax.imshow(fault.dataset[quantity].values, **p)
            ims.append(im)

            if abs(zmin) > abs(maxDepth):
                ax.set_ylim([-abs(maxDepth), zmax])

            ax.set_ylabel('Depth [km]')
            ax.set_xlabel('Time steps')

            # if per‐axis scaling (no global Norm), attach vertical colorbar
            if not use_norm:
                fig.colorbar(im, ax=ax, orientation='vertical',
                             fraction=0.02, pad=0.02)

        # if a global Norm (or explicit vmin/vmax override) was provided, add one shared horizontal bar
        if use_norm and ims:
            cbar = fig.colorbar(ims[0], ax=axes,
                                orientation='horizontal',
                                fraction=0.05, pad=0.1)
            cbar.set_label(quantity)

        return fig, axes


    def PlotSlipRateNew(self,maxDepth=20):


        cmap_coseismic=cmap_tandem.ReturnCmap(vmin=1e-19, vmax=1,Vths=1e-3)
        maxTime=len(self.faults[0].dataset['Time'].values)



        heights=[]
        for i,fault_i in enumerate(self.faults):

            zmin=np.min(fault_i.dataset['z'].values);zmax=np.max(fault_i.dataset['z'].values);
            if np.abs(zmin)>np.abs(maxDepth):
                zmin=-np.abs(maxDepth)
            y_range = np.abs(zmax - zmin)

            heights.append(y_range)


        fig = plt.figure(figsize=(12, sum(heights)))  # Adjust height proportionally to sum(heights)
        gs = GridSpec(len(heights), 1, height_ratios=heights)

        ax=[]
        for i,fault_i in enumerate(self.faults):
            zmin=np.min(fault_i.dataset['z'].values);zmax=np.max(fault_i.dataset['z'].values);y_range = np.abs(zmax - zmin)
            extent=[0,maxTime,zmin,zmax]
            ax_i = fig.add_subplot(gs[i, 0])
            ax.append(ax_i)

            plottingProperties={'vmin':-19,'vmax':0,'aspect':'auto','extent':extent,'cmap':cmap_coseismic}
            dataToPlot=np.log10(np.abs(fault_i.dataset['slip-rate0'].values))

            ax_i.imshow(dataToPlot,**plottingProperties)
            ax_i.set_ylabel('Depth[km]')
            ax_i.set_xlabel('Time steps')
            #fault_i.PlotDataset(ax=ax_i,quantity='slip-rate0',plottingProperties=plottingProperties)
            if np.abs(zmin)>np.abs(maxDepth):
                ax_i.set_ylim([-np.abs(maxDepth),zmax])

        return fig,ax

    def PlotCatalogs(self,ax=None):

        if ax is None:
            fig,ax=plt.subplots()
        else:
            ymin, ymax = ax.get_ylim()


        for i,fault_i in enumerate(self.faults):
            plottingProperties={'color':self.colors[i],}
            fault_i.PlotEventsCatalog(ax=ax,plottingProperties=plottingProperties)


    def PlotFaultPosition(self,ax=None,plottingProp=None):
        if ax is None:
            fig,ax=plt.subplots()

        if plottingProp is None:
            plottingProp={}


        for i,fault_i in enumerate(self.faults):
            ax.plot(fault_i.dataset['x'],fault_i.dataset['z'],label=str(i),color=self.colors[i],**plottingProp)

        #ax.legend()










#%%
def LoadDomain(path,colors=None,file_pattern="fault_*",):
    fault_files=sorted(glob.glob(path+file_pattern))

    faults=[]

    if len(fault_files)==0:
        raise ValueError("couldn't find any files")
    for file_i in fault_files:
        data=xr.open_dataset(file_i)
        faults.append(faultObject(data))

    return domain(faults,colors=colors)
