Fluorescence ray tracing simulation for laser cooling of solids
===
These simulation codes are to simulate the mean fluorescence wavelength and fluorescence escape efficiency of a bulk crystal by by Monte Carlo method.

## Description
- All the parameters for simulation are written in "params.yaml".
- The simulation requires spectroscopic data of the crystal. The spectroscpic data need to be in the directory "spectra" in CSV format. The first column is wavelength in nm, the second column is the normalized fluorescence intensity, and the third column is the absorption cross section in cm^2.
- The file name of the spectroscopic data needs to be named as "{crystal name}\_{temperature}K.csv" for isotropic crystals. For uniaxial crystals, two files for two polarizations are needed and they should be named as "{crystal name}\_{"sigma" or "pi"}\_{temperature}K.csv".


## Requirement
Julia of v.1.8 or later is recommended (we developped the codes under v.1.8.3). For installation of Julia, please visit the homepage of Julia: https://julialang.org/downloads/


## Quick start guide
1. Edit "params.yaml"
2. Run Julia (in shell)
   ```shell
   >julia
   ```
   The following is to run Julia with three additional worker processes, so there will be 4 processes including 1 master process. The code will split the computational task equally to the existing worker processes.
   ```shell
   >julia -p 3
   ```
3. Import packages: CSV, DataFrames, YAML, Distributions, Interpolations
   ```julia
   pkg>add [package name]
   ```
   or
   ```
   julia>import Pkg
   julia>Pkg.add("[packagename]")
   ```
4. Include "main.jl"
   ```julia
   include("main.jl")
   ```
5. Run function main()
   ```julia
   main()
   ```
6. Once the simulation finishes it creates a directory "result" and sub-folder named with the time-stamp containing a CSV file and param.yaml.
7. You can generate a report in pdf by running "report.ipynb".

## Licence
[MIT](https://github.com/hirokitanaka-ikz/fluorescence_ray_tracing/main/LICENCE)


## Reference
[1] Hiroki Tanaka, Stefan Püschel, ...


[hirokitanaka-ikz](https://github.com/hirokitanaka-ikz)
