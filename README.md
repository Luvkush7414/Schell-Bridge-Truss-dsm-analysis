# 2D Pin-Jointed Truss Analysis: Schell Memorial Bridge

This repository contains a structural engineering project that utilizes **Matrix Structural Analysis (MSA) to evaluate an idealized 2D pin-jointed truss model of the historic Northfield Schell Bridge. The entire analysis is programmed from scratch in **MATLAB** using the **Direct Stiffness Method (DSM)**.

## Project Overview

The objective of this project is to programmatically determine the structural response of a complex, multi-span bridge under applied loads. By circumventing commercial GUI-based software and writing the DSM algorithm directly, this project demonstrates a fundamental understanding of linear algebra applications in structural mechanics.

### Key Computational Features
* **Global Stiffness Matrix Assembly:** Programmatic generation of local stiffness matrices $(k')$ and transformation to global coordinates (K).
* **Automated Load Distribution:** Calculation of member self-weights distributed as lumped point loads.
* **Displacement & Force Extraction:** Solving the equilibrium equation $Q = KD$ to extract joint displacements, internal member axial forces (tension/compression), and support reactions.

## Structural Model Definition
* **Total Joints:** 34
* **Total Members:** 75
* **Boundary Conditions:** Pin supports at Joints 1, 17, 33, and 34 (restraining X and Y translation).
* **Material Properties:** * Modulus of Elasticity (E): 200 GPa
  * Cross-Sectional Area (A): 0.01 m²
  * Material Density (rho): 7850 kg/m³

## Applied Loading
1. **Dead Load (Self-Weight):** Computed using W = \rho \cdot A \cdot L \cdot g and applied to connecting nodes.
2. **Dead Load (Deck):** 10 kN downward vertical load applied across bottom chord joints.
3. **Live Load:** 50 kN downward vertical load applied to specific top chord joints (24-28) to simulate vehicular loading.

## Repository Structure
* `/src` - Contains the main MATLAB script (`truss_analysis_main.m`).
* `/data` - Contains the raw joint coordinates and member connectivity data.
* `/docs` - Contains the comprehensive MSA Project Report detailing the assumptions, historical context, and mathematical methodology.

## How to Run the Analysis
1. Clone this repository to your local machine.
2. Open `src/truss_analysis_main.m` in MATLAB.
3. Run the script. The command window will output the resulting member axial forces and support reactions, and a scaled deformation plot will be generated comparing the unloaded and loaded states.
