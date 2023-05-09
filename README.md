Fluorescence ray tracing simulation for laser cooling of solids
===
These simulation codes are to simulate the mean fluorescence wavelength and fluorescence escape efficiency of a bulk crystal by by Monte Carlo method.

## Description
- All the parameters for simulation are written in "params.yaml".
- The simulation requires spectroscopic data of the crystal. The spectroscpic data need to be in the directory "spectra" in the tab-separeted txt. The first column is wavelength in nm, the second column is the fluorescence intensity normalized to 1, the fourth column is the absorption cross section in cm^2.
- The file name of the spectroscopic data needs to be named as "<crystal_name>_<temperature>K.txt" for isotropic crystals. For uniaxial crystals, two files for two polarizations are needed and they should be named as "<crystal_name>_<"sigma" or "pi">_<T>K.txt".


## Requirement
Julia of v.1.8 or later is required. For installation of Julia, please visit the homepage of Julia: https://julialang.org/downloads/


## Quick start guide
1. Edit "params.yaml"
2. Run Julia
   ```shell
   >julia
   ```
   For running Julia in multi-processes. The following is in the case of 3 worker processes, sot there will be 4 processes including 1 master process.
   ```shell
   >julia -p 3
   ```
3. import packages: CSV, DataFrames, YAML, Distributions, Interpolation
   ```julia
   pkg>add [package name]
   ```
4. include "main.jl"
   ```julia
   include("main.jl")
   ```
5. Run function main()
   ```julia
   main()
   ```
6. Once the simulation finishes it creates a directory "result" and sub-folder named with the time-stamp containing a CSV file and param.yaml.
7. You can generate a report in pdf by running "report.ipynb".