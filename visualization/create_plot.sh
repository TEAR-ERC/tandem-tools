# Example code to obtain the plot
# e.g., save_dir=/your/path/to/simulation/result
save_dir='.'

# Example 1. Extract and save essential part from outputs, and returns cumulative slip vs. depth plot, event analyzing plot and image of slip rate in timesteps
python get_plots.py $save_dir -c -wf 24 -csl -dtcr 2 -dtco 0.5 -imsr -ts -evan

# Example 2. Returns cumulative slip vs. depth plot with spin-up after 2.5 m of slip and depth histogram on the side
# python get_plots.py $save_dir -csl -dtcr 2 -dtco 0.5 -spup 2.5 -dd

# Example 3. Returns 4 figures: slip vs. time plot at 7.5 km depth, initial stress profile, and images of slip rate and shear stress in timesteps
# python get_plots.py $save_dir -sl 7.5 -ist -imst -imsr -ts 

