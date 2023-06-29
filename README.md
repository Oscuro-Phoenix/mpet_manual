# main.py

## `run_simulation(config,outdir)` 

This is the key function that takes in a configuration file (`config`) and an output directory (`outdir`) to run the simulation. It heavily relies on all the other files (particularly `mod_cell.py` and `sim.py`) in the codebase to setup the multi-scale system of equations that need to be solved.

1. Get the timescale from config
2. Declare variable for logging data as `log`, declare the solver as `daesolver` , create a simulation object that will be used to initialize and run the simulation `simulation`  and declare a `datareporter`
3. Assign the solver as an LU decomposition sparse solver that is imported from pySuperLU
4. Enable data reporting for all variables except particle and cell ports -- this reporting is set using the datareporting function on `simulation.m` where m refers to the model of the cell that is defined in `mod_cell.py` and is a daeModel (see Line 29 in mod_cell.py).
5. Initialise simulation using solver, reporter and logger using `simulation.initializer()` with appropriate arguments
6. Solve for the DAE variables at t = 0
7. Run the simulation 

## `main(paramfile)`   

This is the main function that set the start time, parses the paramfile into a config object made using `Config()` that is a function in configuration.py in the config folder. Then, the output directory is verified (that it exists) and the paramfile is copied into a file named `input_params_system.cfg`. Also by looping over all the electrodes in the system new `.cfg` files for each trode is populated with information from the paramfile through `config` object. Finally, `config.write()` creates another file that has all the parametric info about the system (look into write function in configuration.py).

Next, the main function checks if there is any git commit info, prints the mpet version, prints git commit info if it was found to exist and this "printing" takes place into a file named `run_info.txt` through the object `fo`. In case the simulation is not a simple galvanostatic simulation with segments of constant current/voltage, the DAETools evaluation tree approach is activated. 

Finally, using the `run_simulation()` function the simulaiton is run, The total time taken is printed.  `shutil` is a python library for copying and archiving files and can be used to move the simulation output around.

# sim.py

## `class SimMPET(dae.daeSimulation)`

The entire file is dedicated to creating this class. Below is a summary of each function 

- `_init_()` - Declares local config object, timescsale variables, by default sets the previouscurrent and applied potential as 0, sets the tolerances and fianlly defines the variable `self.m` that stores all the information about the system of DAE that needs to be solved. `self.m` is set using the  `mod_cell.ModCell()` function that is in a separate file named mod_cell.py that will be discussed later.
- `SetUpParametersAndDomains()`- Uses the information in config to create arrays for each part of the cell namely separator, a (anode) and c (cathode). Within  a and c, additional arrays are created for each particle type (say, particles of different sizes are of different types) (DmnPart) and the number of each particle (for each trode tr we have ith Nvol and jth Particle type for which a fixed number of them are there and this number is based on a particle size distribution `psd_num`).
- `SetUpVariables()` - We will only consider the case where there is no previous data (i.e. `config["prevDir"] = "false"`).  
  - For each trode the avg. filling fraction (`ffrac`) is set to an initial guess. Following this, the volumetric reaction rate R_Vp for each control volume (Nvol) is set to zero as initial guess.  Potential in the solid phase (`phi_bulk`) is set to reference voltage for anode and a voltage for cathode that comes from the config file. Another for loop over j the number of particles in each trode sets the initial concentration as constant (`cs0`) for each particle for type j in control volume i.
  - The potential applied to the cell is initialized based on the type of run (for example if it is a constant voltage run then it comes from `config['Vset']`. If it is not a single volume simulation (`SVsim`) , additional ghost points for discretization in concentration and potential in electrolyte (lyte) are needed and are accordingly set.
  - Finally, the concentration and potential at all other points in the lyte are initialized with an initial guess based on the config file. Additionally, the concentration of electrolyte that interacts with each particle is also set to an intial value `config["c0"]`.  
- `Run(self)` -  There is already a run function that came with the fact that `SimMPET` inherits `dae.daeSImulation` but this function is overloaded with this function that prints simulation progress, integrates the DAE system and keeps checking if an end conditions is reached to stop the simulation. which happens when any of the `self.m.endCondition.npyValues` are nonzero.

# mod_cell.py

## `class ModCell(dae.daeModel)` 

This file is again dedicated to the creation of the ModCell class that encapsulates all the information about the system of equations that needs to be solved. Below is a summary of each function

- `_init_()`  - 
  - First the intiialization function within DAE Tools is run. This is then followed by defining local config, profileType, Number of control vols (Nvol) and Numer of particle types (Npart) and the list of trodes in the system.
  - DmnCell and DmnPart dictionaries are next setup to contain daeDomains for each part of the cell. This is exactly what is later populated with arrays in sim.py.
  - Key simulation variables are defined namely `c_lyte, phi_lyte, phi_bulk, phi_part, Rvp` and `ffrac` . The `daeVariable` function is used to appropriately define the type of variable (either conc_t or elec_pot_t -- see daeVariableTypes.py for more info) and the domain where the variable lives. 
  - Next, if the simulation is not a single variable sim (`SVsim`) then additional ghost point variables are defined. These are again set in sim.py as discussed earlier.
  - Additional variables that are macroscale like `phi_applied`,`phi_cell`, `current` and `endCondition` are defined again using the `daeVariable` function.
  - Ports are setup to help the particles "talk" to the lyte and bulk phases. Essentially, portsOutlyte represent outlet ports that go from lye to trode particles. Similarly, portsOutBulk represent outlet ports that go from bulk to trode particles. This section starts off by setting up empty numpy arrays for the ports which are then populated using the ports module that is imported in the beginning of the code. The particle model is also set here as either being 1variable or 2variable type particle (particle models are setup in mod_electrodes.py). Finally, the ports are connected to the particles (using `ConnectPorts()` function that is inherited from `dae.daeModel`) which already have daeInlet ports that were defined within mod_electrode.py. All of this is being done to ensure that the particle model has locally access to the lye concentration, lyte potential and bulk potential.
- `DeclareEquations()` - 
  - The comment sin code are quite explanatory here but essentially a very standard format is used to declare equations by looping over the respective domain/region and using the `CreateEquation()` function to create an equation for the residual (that has to be zero) using the daeVariables defined earlier.
  - Aside from the standard MPET equations, some peculiar equations are also present for definite the output port variables namely `self.c_lyte, self.phi_lyte` using the portsOutlyte and `self.phi_part` using portsOutBulk. Note: These port variables were defined earlier in `_init_()` (see 3rd sub-bullet)



# mod_electrodes.py

This file is dedicated to writing the equations for the electrode particles which can either be 1 var or 2 var type for each volume element with index vInd and particle index pInd.

## class Mod2var(dae.daeModel)

- `_init_()` - 
  - First the init function from the dae.daeModel is called using super() then self.config, self.trode, self.ind which are config file, trode ID and tuple of volume and particle indices respectively are created. Domain withint he particle is declared as a dae discretization domain where the equations will be defined later. 
  - Next, key variables are defined that include c1,c2,avg. c, avg. c1, avg. c2 and time derivative of avg. c. Additional variables for reaction namely Rxn1 and Rxn2 are defined for each region in particle. Finally, the reactiont type (rxnType) is read off from the config file and a function that can return the reaction rate is defined based on the rxnType this function is `self.calc_rxn_rate()`. This definition happens using the utils module that has the `import_function()` function.
  - Finally, ports that were referred to in last bullet of the `_init_()` function in mod_cell.py are declared so that the particle model has access to the potential in lyte, conc. in lyte and potential in bulk. Additional variables including potential in electrolyte, conc. in electrolyte and potential solid are accordingly declared using the appropriate ports. 
- `get_trode_param()` -  This is a utility function that helps get value of params in config file and will be used in the when the equations are defined in the next function
- `DeclareEquations()` - 
  - Starts with extracting the number of grid pts in the particle, the temperature and the discretization volume fraction vector 
  - Next, if noise has to be added (based on config file), nose1 and noise2 are setup as the 1D interpolation of random normal noise of size numnoise (number of noise pts) and N (number of  grid pts). The noises are stored in a tuple.
  - Two key parameters the chemical potential of the oxidized state (mu_O) and activity of lyte (act_lyte) are calculated using the function calc_mu_O function that can be found near the end of the file
  - Next, 2 equations are setup for the defining the avg filling fraction in the particle and the definition of volume average is used to create a for loop for the residuals (eq1 and eq2). Overall avg. filling fraction is also similarly defined by setting eq.Residual.
  - Two new numpy arrays c1 and c2 are created and the values of self.c1 and self.c2 are copied into it. These two arrays will be sent along with mu_O and act_lyte to the function(s) that set the governing equation in the particle phase.
- `sld_dynamics_0D2var()` -
  - Refer to `sld_dynamics_1D2var()` since this function follows a similar albeit, a simpler structure.
- `sld_dynamics_1D2var()` - 
  - Extract number of grid pts (N) and the temperature (T) 
  - Then the mass matrix for the Discrtized PDE system (Mmat: Mmat*dcdt = RHS) where RHS has the discretized flux terms and the reaction term is obtained using the get_Mmat function around  the end of the file. Also dicretization size dr and edges are obtained from the geometry module (geo).
  - Depending on the type of rxn model we might want to calculate the chemical potential only on the surface or at all grid pts. If it is only on surface (for diffn2, CHR2 models) the surface c1_surf, c2_surf and mu1R_surf, mu2R_surf are used with muO to obtain the overpotentials eta1 and eta2 and subquently obtain Rxn1 and Rxn2 terms.
  - On the other hand, if we have the ACR2 model the self.Rxn1 and self.Rxn2 terms will have to defined on a per grid point basis that is exactly what is done later on the for loop that sets the eq1.Residual and eq2.Residual for each reaction term
  - Next, solid particle fluxes and boundary conditions are populated into the RHS vector for "dffn2" and "CHR2" models. 
  - The final 4 lines create and set the residual for the equations defining the concentration profile in the solid for both phase1 and 2 



# utils.py

This file contains utility functions that are useful throughout the codebase. Linear mean, Harmonic mean, Weighted Linear and Weighted Harmonic mean functions for a vector are first defined. A function for appending ghost points to a vector which was used in mod_cell.py. Other utility functions including get_asc_vec used in mod_cell.py and geometry.py for getting a full numpy array of variables spanning all domains are present. Most of these functions are self explanatory. The last function i.e `import_function(filename,function,...)` is specifically important for the reacrtion rate function, diffusivity function etc. that are material dependent and requires importing the appropriate function from the module that may be present in the mpet package. For example, in the mpet/electrode/materials folder there are several such named functions for chemical potential and activity in of Li in the active material.

# props_am.py

The entire file is dedicated to defining the class `muRfuncs()`. An object of this class essentially has the function that defines the chemical potential of Li in the electrode material. This function is imported using the `import_function()` in utils and is available for use in `calc_muR` function in mod_electrodes.py. Additional material specific helped functions are also defined within the scope of class `muRfuncs()` for example the 1 param graphite model (LiC6_1param.py in electrode/materials) requires `graphite_1param_homog_3()` function while LFP (electrode/materials/LiFePO4.py) requires `reg_sln()` function.

# geometry.py

Has 4 functions namely `get_unit_solid_discr` for getting grid pts and volfrac for the discretized grid in a given shape of solid active particle, `get_dr_edges` to obtain the edge length for area calculation in mod_electrodes.py, `calc_curv` for curvature calculation that is used in certain nonhomgeneous terms in chemical potential of Li in active material and `get_elyte_disc` for linear discretization of the electrolyte phase (this function returns a dictionary). 

# daeVariableTypes.py

Defined the 3 key daeVariable types i.e. concentration, potential and mole fraction.

# data_reporting.py

This is for writing simulation output. For brevity we will only go over `MyMATDataReporter()` class that has a single function `WriteDataToFile` that starts with an empty dict of data that is populated according to `continued_sim` boolean that is set to 1 if we are continuing a previous simulation. Then for each process variable (except for the port variables) a for loop is run where `mdict	` is populated with new data then the `mdict[dkeybase]` value is written into `mat_dat` that lives in the hdf5 file. The only nuance here is that depending on if the `continued_sim` is true additional space needs to be created in `mat_dat` that is done using resize function calls and subsequently assigning the data in `mdict[dkeybased]` to the newly created space. If `continued_sim` is false, fresh space is created using `create_dataset()` function within the object `mat_dat`. 

# Concluding Remarks

`**This section is still in development**`
