Updated on 2024/5/11 Hiroki Ikawa

This is a snowmelt model based on the bulk-surface energy balance approach utilized by Ikawa et al (Water Resources Research, 2024).
bulksurface_sample.jl reads an example data input (usprr2013.csv) created from the data uploaded on Ameriflux, 
and performs model calculations coded in bulksurface.jl.

!! ABOUT model input data
The input file was generated, using usprr_dataset_preprocess.R. There are two important features: 
(1) interpolation of small gaps
(2) calculation of snow precipitation on the floor (psnow) based on snow depth data
