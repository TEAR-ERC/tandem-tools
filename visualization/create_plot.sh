# Example scripts 
# Author: Jeena Yun

# Example 1. Spatiotemporal evolution of sliprate, in timesteps (sliprate_image_timesteps.png)
python get_plots.py examples -fp_prefix outputs/fltst_dp -ev sliprate -step
# Example 2. Spatiotemporal evolution of cumulative slip (cumslip.png)
python get_plots.py examples -fp_prefix outputs/fltst_dp -ev cumslip
# Example 3. Time series of peak slip rate (peak_sliprate.png)
python get_plots.py examples -fp_prefix outputs/fltst_dp -ts sliprate
# Example 4. Time series of shear and normal stresses at 15 km depth (stresses_at_15000m_depth.png)
python get_plots.py examples -fp_prefix outputs/fltst_dp -ts stress -loc 15
# Example 5. Time series of state variable at 5 km depth (state_at_5000m_depth.png)
python get_plots.py examples -fp_prefix outputs/fltst_dp -ts state -loc 5 
# Example 6. Time series of slip at surface (slip_at_surface.png)
python get_plots.py examples -fp_prefix outputs/fltst_dp -ts slip -loc 0
# Example 7. Mesh with cell edges and boundary conditions (mesh_BC.png)
python get_plots.py examples -mesh bp1_sym.msh