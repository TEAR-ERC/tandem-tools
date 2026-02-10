#!/usr/bin/env python3
'''
Create widget to gather Qaukeworx job directory
By Jeena Yun
Last modification: 2026.01.14.
'''

from pathlib import Path
from ipywidgets import Dropdown, Label

"""
Create a dropdown widget to select different plotting options

Returns
-------
dropdown : ipywidgets.Dropdown
    The dropdown widget.
info_label : ipywidgets.Label
    A label that shows the selected option.
get_selected : callable
    Function with no arguments that returns the currently selected option
    (or None if nothing selected yet).
"""

def select_output_dir(**kwargs):
    """
    Create a dropdown widget to select a subdirectory of `base_dir`.
    """
    import os

    if 'base_dir_str' in kwargs:
        # Base directory whose subfolders you want to list
        base_dir = Path(kwargs['base_dir_str'])  # <- change this
    else:
        # If not defined, use default base directory
        base_dir = Path(".")
    default_job_name = 'example_outputs'

    # Gather job names
    job_names = sorted(
        [d.name for d in base_dir.iterdir() if d.is_dir()]
    )
    job_names.insert(0, default_job_name)

    # Do not show pycache for cleaness
    if '__pycache__' in job_names:
        job_names.remove('__pycache__')

    # Create a dropdown menu
    dropdown = Dropdown(
        options=job_names,
        description='Folder:',
        layout={'width': '300px'},
        style={'description_width': 'initial'}
    )

    state = {
        "full_path": base_dir / default_job_name,
        "job_name" : default_job_name,
    }

    info_label = Label(f"Selected job name: {state['job_name']}")

    def on_change(change):
        if change["name"] == "value" and change["new"] is not None:
            state["full_path"] = base_dir / change["new"]
            state["job_name"] = change["new"]
            if state["job_name"] == default_job_name:
                state["full_path"] = '.' / change["new"]
            info_label.value = f"Selected job name: {state['job_name']}"

    dropdown.observe(on_change, names="value")

    def get_selected_job():
        """Return the current selected Path (or None if not selected)."""
        return state['job_name'], state["full_path"]

    return dropdown, info_label, get_selected_job


from ipywidgets import Dropdown, Label

def select_plot_type():
    plot_types = ['', 'evolution', 'timeseries', 'mesh']

    # Create a dropdown menu
    dropdown = Dropdown(
        options=plot_types,
        description='Plot type:',
        layout={'width': '300px'},
        style={'description_width': 'initial'}
    )
    info_label = Label("No type selected yet")

    state = {
        "plot_type": None
        }

    def on_change(change):
        if change["name"] == "value" and change["new"] is not None:
            state["plot_type"] = change["new"]
            info_label.value = f"Selected plot type: {state['plot_type']}"

    dropdown.observe(on_change, names="value")

    def get_selected_plot_type():
        """Return the current selected plot type (or None if not selected)."""
        return state['plot_type']

    return dropdown, info_label, get_selected_plot_type

def select_evolution_type():
    evolution_types = ['', 'State Variable', 'Cumulative Slip', 'Cumulative Slip Contour', 'Slip Rate', 'Shear Stress', 'Shear Stress Change', 'Normal Stress', 'Normal Stress Change']
    args ={
        'State Variable': 'state',
        'Cumulative Slip' : 'slip',
        'Cumulative Slip Contour' : 'cumslip',
        'Slip Rate' : 'sliprate',
        'Shear Stress' : 'shearT',
        'Shear Stress Change' : 'delshearT',
        'Normal Stress' : 'normalT',
        'Normal Stress Change' : 'delnormalT'
    }

    # Create a dropdown menu
    dropdown = Dropdown(
        options=evolution_types,
        description='Evolution type:',
        layout={'width': '300px'},
        style={'description_width': 'initial'}
    )
    info_label = Label("No type selected yet")

    state = {
        "evolution_type": None
        }

    def on_change(change):
        if change["name"] == "value" and change["new"] is not None:
            state["evolution_type"] = change["new"]
            info_label.value = f"Selected evolution type: {state['evolution_type']}"

    dropdown.observe(on_change, names="value")

    def get_selected_evolution_type():
        """Return the current selected evolution type (or None if not selected)."""
        return state['evolution_type'], args[state['evolution_type']]

    return dropdown, info_label, get_selected_evolution_type

def select_timeseries_type():
    timeseries_types = ['', 'State Variable', 'Cumulative Slip', 'Slip Rate', 'Shear Stress', 'Normal Stress']
    args ={
        'State Variable': 'state',
        'Cumulative Slip' : 'slip',
        'Slip Rate' : 'sliprate',
        'Shear Stress' : 'shearT',
        'Normal Stress' : 'normalT',
    }

    # Create a dropdown menu
    dropdown = Dropdown(
        options=timeseries_types,
        description='Time series type:',
        layout={'width': '300px'},
        style={'description_width': 'initial'}
    )
    info_label = Label("No type selected yet")

    state = {
        "timeseries_types": None
        }

    def on_change(change):
        if change["name"] == "value" and change["new"] is not None:
            state["timeseries_types"] = change["new"]
            info_label.value = f"Selected time series type: {state['timeseries_types']}"

    dropdown.observe(on_change, names="value")

    def get_selected_timeseries_types():
        """Return the current selected time series type (or None if not selected)."""
        return state['timeseries_types'], args[state['timeseries_types']]

    return dropdown, info_label, get_selected_timeseries_types

def select_mesh(save_dir):
    from glob import glob

    mesh_path = glob(str(save_dir / '*.msh'))[0]

    return mesh_path

def select_lua_scenario(save_dir):
    from glob import glob

    toml_path = glob(str(save_dir / '*.toml'))[0]

    with open(toml_path, 'r') as fid:
        for line in fid:
            if 'lib' in line.strip():
                lua_lib_name = line.strip().split('\"')[1]
            if 'scenario' in line.strip():
                lua_scenario_name = line.strip().split('\"')[1]

    return lua_lib_name, lua_scenario_name