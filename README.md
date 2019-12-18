
# A simple river-routing model written in MATLAB

It takes 0.5 degree run-off data from 5 climate models and distributes them over a 6-minute resolution Hydro STN netwrok.
The model routs the runoff over the river netwrok to create a simple river flow dataset.

Simplifications: 
- Time step is 1 day (equal to resolution of the run-off data)
- Slope of network is assumed to be even across the glob; each grid-cell flows into the next gridcell every time step, regardless of the slope

* GCM run-off data are not uploaded due to Github's upload limitations
* The produced river network plots are in the Figures directory

