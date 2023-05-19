# Example code to obtain the plot
save_dir=/your/path/to/simulation/result

# Example 1. Returns cumulative slip vs. depth plot with spin-up after 2.5 m of slip and depth histogram on the side
python get_plots.py $save_dir -c -csl -dtcr 2 -dtco 0.5 -spup 2.5 -dd

# Example 2. Returns 4 figures: slip vs. time plot at 7.5 km depth, initial stress profile, and images of slip rate and shear stress in timesteps
python get_plots.py $save_dir -sl 7.5 -ist -imst -imsr -ts 

