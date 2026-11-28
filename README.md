# Vibration Analysis of an Agricultural Tractor (2-DOF Model)

A numerical and analytical study of the vertical bounce and pitch vibration dynamics of an agricultural wheeled tractor subjected to road excitation. Developed as part of **ME3103: Dynamics of Machinery** (Team 5), this repository contains the complete MATLAB simulation suite, frequency-domain Power Spectral Density (PSD) analysis, phase-plane stability portraits, kinematic animations, and the final project presentation.

---

## Table of Contents

1. [Repository Structure](#repository-structure)
2. [Mathematical Formulation](#mathematical-formulation)
   - [Governing Equations of Motion](#governing-equations-of-motion)
   - [System Parameters](#system-parameters)
   - [Free Vibration Analysis](#free-vibration-analysis)
3. [Simulation Architecture](#simulation-architecture)
4. [Results & Visualizations](#results--visualizations)
   - [Frequency Sweep (Chirp) Response](#1-frequency-sweep-chirp-response)
   - [Power Spectral Density (PSD) Analysis](#2-power-spectral-density-psd-analysis)
   - [Phase-Plane Stability Portraits](#3-phase-plane-stability-portraits)
   - [Real-Time Kinematic Animation](#4-real-time-kinematic-animation)
5. [How to Run](#how-to-run)
6. [References](#references)

---

## Repository Structure

```text
.
├── FINAL_PPT_ME3103_TEAM5.pptx            # Final project presentation slides
└── MATLAB_stuff/
    ├── code/
    │   ├── tractor_simulation_chirp.m     # Swept-sine (0–10 Hz) time-domain simulation & RMS metrics
    │   ├── tractor_psd_analysis.m         # FFT-based Power Spectral Density (PSD) analysis
    │   ├── tractor_phase_portraits.m      # Vertical & pitch phase-plane trajectory generator
    │   └── tractor_simulation_kinematics.m # 2-DOF tractor kinematic animation generator
    └── results/
        ├── tractor.png                    # Tractor reference / model illustration
        ├── tractor_flowchart.png          # Numerical integration & solver architecture flowchart
        ├── CHIRP_graphs.png               # Vertical bounce & pitch response (0–10 Hz sweep)
        ├── SWEPTSINE_results.png          # Command window RMS acceleration outputs
        ├── PSD_graphs.png                 # Vertical & pitch Power Spectral Density plots
        ├── Tractor_Phases.gif             # Animated phase-plane stability trajectories
        └── Tractor_Kinematics.gif         # Animated 2-DOF tractor chassis motion
```
$$
\mathbf{M}\ddot{\mathbf{q}} + \mathbf{C}\dot{\mathbf{q}} + \mathbf{K}\mathbf{q} = \mathbf{F}
$$
