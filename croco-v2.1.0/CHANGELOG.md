# Changelog

Release changelog are available here : https://gitlab.inria.fr/croco-ocean/croco/-/releases

## [2.1.0] - 2025-03-21

### Added

- NEW TEST FRAMEWORK : add benchmarking python tools
  to perform tests runs, see documentation or BENCH/README.md for more details.
  Note that cppkeys containing RVTK have been renamed using CVTK (Issue #298)

- NEW TEST CASE : FLASH_RIP (Issues #291 and #323)

- NEW DIAGNOSTIC : add mixed layer diagnostics computation (criterion : 
  density, temperature, bvf) (Issue #218)

- NEW CPPKEY : add RAIN_FLUX cppkey to take into account water and temp flux 
  from rain (Issue #195)

- MUSTANG : add dredging feature (Issue #250)

- COUPLED RUN SCRIPTS

  - Add possibility to give a pre-built grid file to OASIS
  - Add new script to take into account the ideal WRF configurations
  - Add possibility to have different WRF output frequencies between domains

### Fixed

- IO
  - Fix error in wrt_his.F when nrechis=0 and nrpfhis=1 (Issue #172)
  - Fix wrong array name for XIOS and WKB_WAVE (Issue #260)
  - Fix wrong hbl test on XIOS field activation with GLS_MIXING (Issue #264)
  - Fix FILLVAL, mask was misdone on variable value 
    instead of mask value (Issue #324)
  - Fix writing of sediment layers (SEDIMENT) (Issue #329)
  - Fix rstTime netcdf index mismatch between PISCES and CROCO (Issue #331)

- COMPILATION
  - Fix compilation error with TVD and NBQ (Issue #191)
  - Fix NBQ+XIOS compilation (Issue #285)
  - Fix parallel compilation (except for AGRIF) (Issue #309)

- COUPLED RUN SCRIPTS
  - Back to previous default option in create_config 
    interactive q/a concerning matlab tools (Issue #239)

- RIVER
  - Fix PSOURCE_MASS capabilities broken by previous change (Issue #252)
  - Fix MPI reproducibility with PSOURCE when a source is 
    in one MPI domain rejecting towards another MPI domain (Issue #305)

- NUMERIC
  - Fix misuse of temporary WFe,WFx arrays for horizontal w 
    advection in NBQ (Issue #258)
  - Fix typo in order 5 scheme (Issue #318)

- NBQ
  - Fix MPI reproducibility with NBQ (Issue #316)

- ABL1D (Issues #272, #279 and #321)
  - Fix sea surface currents in ABL1D
  - Fix wrong mixing length computation
  - Fix ABL1D perfect restart
  - Fix ECUME6/ABL1D
  - Fix cppdefs organisation and incompatibility with update of ABL1D

- MUSTANG
  - Fix several bugs with MUSTANG output and initialisation (Issue #267)
  - Add a minimim porosity in deposit in V1 (Issue #267)
  - Fix transition with ero_option=3 (Issue #267)
  - Fix MPI reproducibility using MUSTANG_CORFLUX (Issue #319)
  - Fix SUBSTANCE_SUBMASSBALANCE feature on river fluxes (Issue #322)

- OTHERS
  - Fix initialisation of non-declared arrays (Issue #261)
  - Clean unused variables in set_diags_ek.F and set_diags_pv.F (Issue #263)
  - Clean comment in sta.h (Issue #299)
  - Fix and clean compatibility AGRIF and TIDE_MAS (Issue #308)
  - Fix data instruction replace by a declaration 
    of array and initialisation of it (Issue #311 and #328)
  - Fix step3d_fast with wave but no wet and dry (Issue #332)
  - Fix wrongly place mask in u2bc (Issue #338)

### Changed

- COMPILATION
  - Jobcomp script can now be executed with terminal options, see 
  ```./jobcomp -h``` for more details (Issue #298)

- WAVE
  - In MRL_WCI change limits of wave height in surfzone (Issue #306)
  - Change wavemaker (periodic boundaries + single-sum), 
    check of mean angle after the procedure and periodization is 
    abandoned if the mean angle change is too high (Issues #291 and #323)

- MUSTANG
  - Change for MUSTANG output (Issue #163)
    - Avoid possibility of overlapping in vname by 
      using a separate array vname_must
    - Use l_out_subs from substance namelist to allow the output of only wanted 
      substance
    - Add boolean in namelist for choosing which variables to output and 
      allocate only the needed arrays
    - Update paraMUSTANG_defaults.txt with new booleans availables
    - Update XIOS output to have the save available variables in all MUSTANG
      output options
  - Change for MUSTANG initialisation, review of namelist parameters (Issue #250)

- PISCES 
  - Improve and add changes (Issue #141 and Issue #281)
    - Phase the PISCES version with the one used for the CMIP7 exercise (NEMO 4.2.*)
    - Update on the PISCES interfacing module between NEMO and CROCO
    - Improve Diagenetic module : performance and diagenetic processes (e.g. increased number of POC classes, ...)
    - Add creation of an independent pisces restart file (managed in namelist_pisces_ref) to improve restartability
    - Rename the simplified version of PISCES, cpp key pisces_npzd
    - Fix some bugs
    - Optimize PISCES code on
      - representation of the lability of the particle pool
      - several optimizations to the calculation of certain variables (performance).
  - Reformulate slightly the microzooplankton grazing as a function
    of the prey size (Issue #340)

- COUPLED RUN SCRIPTS
  - Improve possibility to add weight file for OASIS
  - Remove WRF rsl files because too heavy 
  - Modify atm_his_h in hours as atm_his in minutes
  - Update for suporting both WRFV4.2.1 and WRFV4.6 OASIS names
  - Update create_config for including pytools

### Deprecated

/

### Removed

- Remove CPPKEYS :
  - key_MUSTANG_specif_outputs, key_MUSTANG_add_consol_outputs 
    (see Issue #163 MUSTANG outputs are now all specified by namelist)
  - key_CROCO, MORPHODYN_MUSTANG_byHYDRO
  - RVTK, RVTK_DEBUG, RVTK_DEBUG_ADVANCED have been renamed (RVTK changed 
    to CVTK)

- Remove unused file OCEAN/spkitlocal_nh.F90 (Issue #259)

- Remove CVTK previous test framework

### Other

- Add information about the branch status (stable or not)
  in the README file (Issue #271)

### Contributors on this release

- Contributors already on board : 
  Rachid Benshila, Maurice Brémond, Matthieu Caillaud, Gildas Cambon, 
  Nicolas Ducousso, François Dufois, Christian Ethé, Jonathan Gula, 
  Swen Jullien, Solène Le Gac, Florian Lemarié, 
  Patrick Marchesiello, Camille Mazoyer, 
  Cyril Nguyen, Renaud Person, Joris Pianezze

- New contributors : 
  Lisa Maillard, Anne-Lou Schaefer, Simon Treillou, Sebastien Valat