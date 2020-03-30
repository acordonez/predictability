# predictability
Scripts for the paper [Processes Controlling Arctic and Antarctic Sea Ice Predictability in the Community Earth System Model](https://journals.ametsoc.org/doi/full/10.1175/JCLI-D-18-0348.1)

These scripts contain the analysis Matlab code used as the basis of this paper.

The model runs used are control runs from the CESM large ensemble project. Note that we do not use the 
historical period ensemble, just the 1850's controls. More information about these simulations
is available here: http://www.cesm.ucar.edu/projects/community-projects/LENS/. The script 
'lag_correlations_LENS_modern.m' will generate correlations from the historical runs, though 
these are only briefly discussed in our paper.

These scripts were written and run on the NCAR Yellowstone and Cheyenne computers / Geyser cluster.

For questions and comments, please contact Ana Ordonez at anaordonez@comcast.net

## Scripts for Figures
Fig. 1: ice_area_hemi_seasonal_plot_fig1.m  
Fig. 2: ice_extent_regional_fig2.m  
Fig. 3: not available; similar to Fig. 4
Fig. 4: monthly_lag_correlation_plots_fig4_5.m  
Fig. 5: monthly_lag_correlation_plots_fig4_5.m  
Fig. 6: September_lead_correlations_CESM_fig6.m, Sept_mar_line_correlation_plots_CESM_fig6.m  
Fig. 7: September_lead_correlations_DYN_regional_fig7.m, Sept_mar_line_correlation_plots_DYN_regional_fig7.m  

## Other functions to note
'grid1togrid2.m' is a regridding script for CESM output. 
'map_gx1v6NH_TO_stereo25km_blin.161123.nc' and 'map_gx1v6SH_TO_SHstereo25km_blin.161213.nc' are regridding weight files generated by the SCRIP utility

'get_monthly_sie_from_obs.m' returns the ice extent from satellite observations. The satellite 
observations are not included but can be obtained from https://nsidc.org/data/nsidc-0079.

Certain colormaps and netcdf utilities may not be included but should be obtainable from Mathworks if 
desired (for example, lbmap is the Light Bartlein color map).

Not all of the region map scripts are included. 'sio_2016_mask.nc' is the Sea Ice Outlook region mask.
The other regional masks are defined within scripts.
