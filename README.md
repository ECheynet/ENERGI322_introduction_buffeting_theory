# ENERGI322_introduction_buffeting_theory
Matlab code used in the lecture introducing the buffeting theory

## Summary

This repository contains Matlab code used in the lecture series ENERGI322 introducing the buffeting theory. The buffeting theory is used in wind engineering to study the dynamic wind-induced responses of structures. The materials provided here aim to offer practical insights into applying the theory using computational methods.

## Content

The repository is organised into two main components:

- **Documentation.mlx**: A Matlab Live Script that serves as the primary documentation. 

- **Matlab Functions Folder**: This directory contains essential Matlab functions used throughout the lecture. The functions include:
  - `randomProcess`: Generates turbulent wind velocities.
  - `getResponseFD`: Calculates the frequency domain response of a structure subjected to wind loading, incorporating aerodynamic and mechanical admittance functions.
  - `getResponse1D_ODE45`: Solves the differential equations representing the motion of a Single Degree of Freedom (SDOF) system under wind loading in the time domain.
  - `applyAAF_TD`: Applies the Aerodynamic Admittance Function (AAF) to wind velocity time series data, adjusting the data to account for the frequency-dependent aerodynamic effects.

The provided code and documentation are intended for educational purposes. Users are encouraged to modify and extend the code to suit their specific needs. This is the first version of the code, some bugs may still be present.
