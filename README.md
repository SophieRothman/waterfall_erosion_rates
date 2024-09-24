MATLAB Codes provided to support “Waterfalls alter reach-scale fluvial erosion rates: Evidence from field data and process modeling” by Sophie D Rothman, Joel S Scheingross and Scott W McCoy, currently in review at JGR: Earth Surface.
[width_calc_tribs_forsup.m] = creates estimates of channel width given lidar data, using the TopoToolbox2 toolkit. Widths from different elevations should be used to calculate a relationship between width and depth for the reach of interest. This relationship is neded to use the ‘combined_processes_model_wrapper_forsupplement.m’ code to estimate combined processes model predictions of reach-scale erosion rates. 
[combined_processes_model_wrapper_forsupplement.m]=Combined processes model calculator - this code is a wrapper that takes inputs from other codes and makes a time averaged combined processes model prediction as well as a planar channel prediction.

The following are codes are called by [combined_processes_model_wrapper_forsupplement.m]:
[all_erosion_sites.m]=Calculates channel flow conditions for a waterfall-rich channel in a given flood level and uses those results to run  plunge pool and total load erosion codes. 
[only_tl_sites.m]= Calculates channel flow conditions for a waterfall-free channel in a given flood level and uses those results to run the  total load erosion code.
[Qsc_pool_public_tw.m]=Calculates sediment flux through waterfall plunge pools, following Scheingross & Lamb, 2016.
[pool_erosion_for_sophie_tw.m]=Uses sediment flux through plunge pools to calculate erosion rates in plunge pools following Scheingross & Lamb, 2017..
[total_load_code_fct_public.m]=Calculates erosion rates in a planar channel due to abrasion following Lamb et al., 2008. Written by Michael P. Lamb and editted by Joel S. Scheingross. 
[c_from_Qs_fct_tw.m]=Calculates the concentration of sediment around plunge pool boundaries following Scheingross & Lamb, 2016.

Additionally, the excell spreadsheet contains the supporting dataset for “Waterfalls alter reach-scale fluvial erosion rates: Evidence from field data and process modeling”, including waterfall, reach-scale and basin-averaged metrics, erosion rates, and pebble count data.
