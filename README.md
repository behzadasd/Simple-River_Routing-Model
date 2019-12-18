
A simple river-routing model written in MATLAB

It takes 0.5 degree run-off data from 5 climate models and distributes it over a 6-minutes resolution Hydro STN netwrok.
The model routs the runoff over the river netwrok to create a simple river flow dataset.

Simplifications: 
- Time steps are 1 day 
- Slope of network is assumed to be even across the glob; each grid-cell flows into the next gridcell every time step, regardless of the slope


