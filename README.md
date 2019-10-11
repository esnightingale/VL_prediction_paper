Author:         Emily S Nightingale
Institutions:   London Schoool of Hygiene and Tropical Medicine, London, UK
Date Published: October 2019

Project:        Spatio-temporal Approach to Short-Term Prediction of Visceral Leishmaniasis Diagnoses in India

This directory contains an analysis of monthly, block-level VL case counts across the endemic region of Bihar and Jharkhand states in North-Eastern India.

Prerequisites
The analysis was done in RStudio Version 3.5.1 - https://www.rstudio.com/. The following packages are required: 

surveillance, hhh4addon, dplyr, plyr, lubridate, data.table, ggplot2, reshape2, spdep, rgdal, maptools, rlist, 
lattice, MASS, pracma, stringr, png, abind, ggrepel, cowplot, gridExtra, fanplot.

Dataset 
Data are owned by the Indian government via the National Vector Borne Disease Control Programme and we do not 
have permission to share them publicly.
 
File list
The following files 
	1. Code scripts used to run the analysis and generate results
		a) 1_data_run.R - Load data and set up R objects for analysis 
		b) 2_functions.R - Run user-defined functions for model assessment, plotting and forecasting
                c) 3_random_model_draw.R - Generate and evaluate random set of models as described in manuscript
		c) 4_selection.R - Systematic selection process (by ranked probability score) over possible model structures within the                                    surveillance framework
		d) 5_evaluation.R - Calculate alternative fit and prediction metrics for all models fit during selection
		e) 6_figures - Produce figures included in manuscript.
		f) plot_hhh4_amended.R - A variation of the automatic plot for hhh4 objects which is included in surveillance

