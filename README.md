# Approximate Bayesian Input Signal Design (MATLAB)

This repository contains MATLAB implementations of the algorithms developed in the paper
> **An Approximate Bayesian Approach to Optimal Input Signal Design for System Identification**  
> Piotr Bania, Anna Wójcik, AGH University of Krakow (2025)
> 
> doi: 10.20944/preprints202509.1390.v2

## 📖 Overview
The design of informative input signals is crucial for accurate system identification.  
Traditional Fisher-information–based approaches are inherently local and may fail under large parameter uncertainty or nonlinear dynamics.  
This project implements a **Bayesian input design framework** based on maximizing the mutual information (MI) between model parameters and observations, using a **tractable lower bound**.

Key contributions implemented here:
- efficient computation of the Kolchinsky–Tracey MI lower bound
- a highly efficient, Kalman filter-based algorithm that makes the method feasible for long experiments by avoiding the inversion of large covariance matrices 
- example scripts replicating selected results from the paper.

## ⚙️ Functions
- `dij_calc.m` - calculation of d<sub>ij</sub> quantity according to Lemma 3 and 4.
- `get_abcg.m` - a function computing discrete-time state-space matrices for an optically pumped magnetometer (section 6.3).
- `get_lfun.m` - calculation of the objective function for finding the MAP estimate. 
- `information_lb.m` - calculation of information lower bound by Lemma 1, 3, 4.
- `information_lb_12.m` - simplified version of the information lower bound for two possible values of the parameter theta.
- `simulate_sys.m` - a quasi-linear system simulator.
