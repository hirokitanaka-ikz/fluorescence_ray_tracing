Monte Carlo fluorescence ray tracing simulation for laser cooling of solids
===
These codes are to simulate the mean fluorescence wavelength and fluorescence escape efficiency of a bulk crystal by the Monte Carlo method. This repository is linked to the following article: Hiroki Tanaka and Stefan Püschel, "Monte Carlo fluorescence ray tracing simulation for laser cooling of solids," Optics Express 32(2), 2306-2320 (2024) Optics Express. https://doi.org/10.1364/OE.503250


## Description
- All the parameters for simulation are written in 'params.yaml'. You can specify the crystal, pump beam, and simulation settings.
- The simulation requires spectroscopic data of the laser-cooling medium. The spectroscopic data need to be in the directory "spectra" in CSV format. The first column is wavelength in nm, the second column is the normalized fluorescence intensity, and the third column is the absorption cross section in cm^2. This repository currently include the spectroscopic data of Yb:LiYF4 (YLF) for a wide range of temperatures (50 - 300 K) and that of Yb:YAG at 300 K.
- The file name of the spectroscopic data needs to be named as '{crystal name}\_{temperature}K.csv' for isotropic crystals. For uniaxial crystals, two files for two polarizations are needed and they should be named as '{crystal name}\_{"sigma" or "pi"}\_{temperature}K.csv'. The codes will be updated for simulating biaxial crystals in the near future.
- This simulation uses absorption cross section values derived from the fluorescence spectra using the reciprocity method, or McCumber theory because the values at long wavelengths are more reliable than the values calculated from the measured transmission spectra using the nominal doping levels.
- Although the shape of absorption cross section spectra derived by the reciprocity method mostly shows a perfect agreement with those measured in transmission spectroscopy, we often observe a discrepancy in their values. Therefore, the correction coefficient 'correction_coeff' in the parameter file 'params.yaml' needs to be adjusted. In the case of Yb:YLF, we found a correction value to be 1.15 (the values from measurement are larger by a factor of 1.15 than the values derived by the reciprocity method).


### Simulation of mean fluorescence wavelength
- If you only want to compute the mean fluorescence wavelength, we recommend you to set the maxium number of reflections to reduce the computational cost ('max_ref_count') to 100 or 1000. Note that, higher the refractive index, larger the required maximum reflection counts to obtain the reliable results.
- In this case, you can put QE=1.0 (internal quantum efficiency to be unity) and alpha_b=0.0 (background absorption coefficient to be zero).
- For a good convergence of the computed values, simulating 10^6 fluorescence rays is recommended.


### Simulation of fluorescence escape efficiency
- If your target is to compute the fluorescence escape efficiency, you should set appropriate values of QE and alpha_b, and the maximum reflection count (max_ref_count) must be infinity (set 0 in 'params.yaml'). Note that the simulation takes long time when there is no limit on max_ref_count.


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
3. Import packages (only if not yet installed): CSV, DataFrames, YAML, Distributions, Interpolations
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
[1] Hiroki Tanaka and Stefan Püschel, "Monte Carlo fluorescence ray tracing simulation for laser cooling of solids," Optics Express (under review).


[hirokitanaka-ikz](https://github.com/hirokitanaka-ikz)
