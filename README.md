# Group Number 12 - ASEN 2004, Section 011 - Bottle Rocket Simulation Code 

## This Github repository houses all the code for the ASEN 2004: Vehicle Design's Rocket Bottle Lab. This code reflects the work of group 12, and should not be copied for usage. This repository provides easy access to view all the files. Short descriptions of the files are given below, along with background information on the project. 

[ Thermodynamic Model - Model 1 ]
- rocketMain_thermo.m : Main driver script for the thermodynamic model, simulates a rocket trajectory with given initial conditions, plots
- state1func.m, state2func.m : specific state functions for ode45 run the actual simulation mathmatics
- opt1.m, opt2.m : stop condition functions to control rocket flight states in ode45

[ Isp Model - Model 2 ]
- rocketMain_isp.m : Main driver script for Isp model, loads static test stand data and finds deltaV for model trajectory
- state1funcROCKET.m : state function for isp model ode45 simulated run, no thrust, instantaneous velocity 

[ Direct Interpolation Model - Model 3 ]
- rocketMain_DI.m : Main driver script for the DI model, combination of both thermodynamic model and ISP
- state1funcDI.m, state2funcDI : state functions based on thermodynamic model with interpolated thrust at ode45 timesteps

[ Monte Carlo Simulations ] 
- monte.m : normally distributed initial conditions, runs through models 1 & 2 with set number of trials
- monteMain.m : runs monte carlo simulation and plots error ellipses, returns averages to estimate rocket launch within percentile categories
- fileLoad.m : loads a file, specific to static test stand data, converts to newtons
- dataRead.m : reads multiple static test stand file data, modifies with tolerence slope check, zero removal, and 4x scaled shift on mean                  negative values, performs limited statistical analysis on fitted isp models

[ Sensitivity Analysis ]
- SensAnalysis.m : performs sensitivity analysis on certain alterable initial conditions, plots to find the best possible initial condition

## Background Information

The purpose of this code was to develop robust simulations that aimed to predict a real life bottle rocket launch. Three distinct models were developed to get a better sense of model limitations. The thermodynamic model used completely theoretical calculations, whereas the Isp model used data from static test stand launches during lab sections. The direct interpolation model combined both of these models to account for benefits of both. Each model had underlying assumptions that deviated it from natural behavior seen during rocket launches. A monte carlo simulation developed from the models allowed the team to provide a solid pre launch estimate for the final distance, crossrange angle, time in the air, and max height. The monte carlo simulation gave confidence intervals of where the rocket would launch, allowing the team to account for normally distributed initial conditions. Sensitivity Analysis allowed the team to highlight a specific design choice that would give the rocket a farther range. Despite the assumptions and simplified behavior, the code estimated the exact final location of the rocket launch within 2 meters in distance, 5 degrees in crossrange angle, and 1 second in time. 
