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

This is the main function that set the start time, parses the paramfile into a config object made using `Config()` that is a function in configuration.py in the config folder. Then, the output directory is verified (that it exists) and the paramfile is copied into a file named `input_params_system.cfg`. Also by looping over all the electrodes in the system new `.cfg` files for each trode is populated with information from the paramfile through `config` object. Finally, `config.write()` creates another file that has all the paramteric info about the system (look into write function in configuration.py).

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
