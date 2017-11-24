# DynamicPhotosynthesis

Source files required to reproduce analysis from the paper "Dynamic modelling of limitations on improving leaf CO2 assimilation under fluctuating irradiance" including the source code of the model and the R scripts required to reproduce all figures in the paper.

Before making use of the scripts, please run the `INSTALL.R` script to install all dependencies. This will require an Internet connection. Also, if you are using Windows, make sure to have [Rtools](https://cran.r-project.org/bin/windows/Rtools/) installed. If you are encounter any issues running this code, please use the Issues tab of this github website to report them.

# Scripts for simulations

The scripts required to reproduce the simulations are located in the folder `Code` organized according to the different types of calculations required. The results of the calculations are cached as binary R files in the folder `Intermediate`, whereas the final figures will be located in the folder `Output`.

# Source code of the model

# Source code of the model

The model was implemented using the [ODEDSL](https://github.com/AleMorales/ODEDSL.jl) library in Julia, and converted into an R package using the [SimulationModels](https://github.com/AleMorales/SimulationModels.jl) R package. This R package is included in the `Packages` folder under the name `ThylakoidMetabolism`. Performing simulations with the model only requires this R package (plus all dependences installed by the `INSTALL.R` file). The actual implementation of the model can be found in the file `Packages/MiniModel/src/MiniModel.cpp`. This file was automatically generated from a high-level description of the model (see below) so no documentation or comments are provided.

The model was implemented using an ad-hoc domain specific language written in the Julia programming language via the [ODEDSL.jl](https://github.com/AleMorales/ODEDSL.jl) library. Thus, although the model may be extended by modifying the source code inside the ThylakoidMetabolism package provided here, it is more convenient to use ODEDSL and generate a new R package with the new version. For that purpose, the `MiniModel.ode` file is provided. The virtual mutants used in Figure 7 of the pap, as well as the modified model used 