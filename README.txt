Author:         Emily S Nightingale
Institutions:   London School of Hygiene and Tropical Medicine, London, UK
Date Published: October 2019

Project:        Spatio-temporal Approach to Short-Term Prediction of Visceral Leishmaniasis Diagnoses in India

This directory contains an analysis of monthly, block-level VL case counts across the endemic region of Bihar and Jharkhand states in north-eastern India.

Prerequisites
The analysis was done in RStudio Version 3.5.1 - https://www.rstudio.com/. The following packages are required: 

surveillance, hhh4addon, dplyr, plyr, lubridate, ggplot2, reshape2, spdep, rgdal, maptools, rlist, 
lattice, MASS, pracma, stringr, png, abind, ggrepel, cowplot, gridExtra, viridis.

Dataset 
Data are owned by the Indian government via the National Vector Borne Disease Control Programme and we do not 
have permission to share them publicly. This repository contains a dummy dataset simulated from the final model fit.
 
File list
The following files are contained in the repositiry:
	1. Code scripts used to run the analysis and generate results
		a) 1_data_run.R - Load data and set up R objects for analysis 
		b) 2_functions.R - Run user-defined functions for model assessment, plotting and forecasting
		c) 3_random_model_draw.R - Fit randomly-constructed models from a subset of realistic model components.
		d) 4_selection.R - Systematic selection process (by ranked probability score) over possible model structures within the surveillance framework
		e) 5_evaluation.R - Calculate alternative fit and prediction metrics for all models fit during selection.
		f) 6_figures.R - Paper figures
		g) plot_hhh4_amended.R - A variation of the automatic plot for hhh4 objects which is included in surveillance. 
	2. Datasets
		a) input_sim.RData - Input dataset simulated from the final model fit.

	3. Variable codebook - describing variables contained in the dataset including data types and value labels
