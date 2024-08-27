# Multidimensional RF Pulse Design with Consideration of Concomitant Field Effects
 
Code for generating simulation results for the original and proposed methods in this paper: *Ziwei Zhao, Nam G. Lee, Krishna S. Nayak. "Multidimensional RF Pulse Design with Consideration of Concomitant Field Effects"*. Submitted to Mag. Reson. Med. 2024

## Set up

### Initialization

#### set up the path
 - Please modify line 3 of `./sim/setup_path.m` to declare the root path where the main folder `multi_RF` is located.

#### compile .c files 
 - Please compile `bloch.c` and `bloch_maxwell.c` scripts under `./third_party/Bloch_simulator` by using `mex` command in the MATLAB commanding window.


## Code Structure
 
### Figures: 
- To reproduce Figure 2, please run `./figures/sim_results/Fig2_validation/mrm_figure2.m`.
  - To generate the excitation profiles (*.mat* files), please use MATLAB script `./sim/sim_concomitantfield_approx.m`
- To reproduce Figure 3, please run `./figures/sim_results/Fig3_distance/mrm_figure3_distance.m`.
  - To generate the excitation profiles (*.mat* files), please use MATLAB script `./figures/sim_results/Fig3_distance/generate_matfiles_scripts_fig3.m`
- To reproduce Figure 4, please run `./figures/sim_results/Fig4_field_strength/mrm_figure4_fieldstrength.m`.
  - To generate the related excitation profiles (*.mat* files), please use MATLAB script `./figures/sim_results/Fig4_field_strength/generate_matfiles_scripts_fig4.m`
- To reproduce Figure 5 (A), please run `./figures/sim_results/Fig5_phantom_055T/mrm_figure5A_tailored_b0b1.m`.
  - To generate the related excitation profiles (*.mat* files), please use MATLAB script  `./figures/sim_results/Fig5_phantom_055T/generate_matfiles_script_fig5.m`
- Supporting Figure S1 was generated inside the function `./sim/sim_concomitantfields_8channel.m` with the setting `display_figure = true`.
- Supporting Figure S2 was generated using `./figures/sim_results/FigS2_duration/mrm_sup_figure2_duration.m`.
  - To generate the related excitation profiles (*.mat* files), please use MATLAB script `./figures/sim_results/FigS2_duration/generate_matfiles_scripts_sup_fig2.m`

### Main simulation functions:
- `sim_concomitantfield_approx.m` generates concomitant field accuracy results (Figure 2).
- `sim_concomitantfields_8channel.m` computes the excitation profiles using original and proposed methods at different isocenter, different main field strengths, different T2 values, and different undersampling foctor of the designed RF pulses using 8 channel setups. The resulting excitation profiles are saved as *.mat* format.
- `sim_concomitantfields_1channel.m` computes the excitation profiles using original and proposed methods at different isocenter using 1 channel setups. The resulting excitation profiles are saved as *.mat* format.

### B0, B1+ and GIRF measurement datasets: 
- GIRF measurements at 0.55T are provided in `./b0b1_map_055T/GIRF_20200221_Duyn_method_coil2.mat`
- B0, B1+ maps with masks at 0.55T with different isocenter (0, 5, 10, 15cm) are provided in `./b0b1_map_055T` folder.
- `./b0b1_map_055T/gene_b0_map_flash3d.mat` is used to estimate a B0 map from flash_2d images using linear fitting method.

### 2D Spin echo sequence:
- `./seq/demo_pulseq_Pauly_1989_JMR_modified.mat` generates the 2D single slice spin echo sequence in Pulseq *.seq* file with the excitation RF replaced with the designed RF pulse.

### Other functions: 
- `./sim/funcs/STA_maxwell_system_matrix_con.m` calculates a small-tip-angle system matrix with the consideration of concomitant fields during the iterative RF design (proposed method).
- `./sim/funcs/STA_maxwell_system_matrix.m` calculates a small-tip-angle system matrix without the consideration of concomitant fields during the iterative RF design (original method).
- `./sim/funcs/figure_out_transformation_matrix.m` and `./sim/funcs/geometry` folder outputs the coordination transformation matrix that is used in the application of GIRF at 0.55T.
- `./sim/funcs/apply_GIRF_tx.m` calculates the GIRF-corrected trajectory.
- `./sim/funcs/calcgradinfo.m` calculates the kspace trajectory from the input gradient waveforms.
- `./seq/input_pulseq_cartesian_spin_echo_multi_RF.m` saves sequence parameters that is used in the pulseq file.
- `./seq/calculate_cartesian_spin_echo_imaging_parameters.m` calculates basic sequence parameters and define gradient performance at 0.55T that are used in the pulseq file.


## Dependencies
[reVERSE_GIRF](https://github.com/mriphysics/reverse-GIRF?tab=readme-ov-file) by Shaihan Malik, July 2015.

[Phase_relaxed_CPMG_excitation](https://github.com/mriphysics/phase_relaxed_CPMG_excitation) by Shaihan Malik, July 2015.

[lsqrSOL](https://github.com/areslp/matlab/tree/master/lsqrSOL) for solving a linear non-square matrix problem. 

[Bloch_simulation](http://mrsrl.stanford.edu/~brian/blochsim/) by Brian Hargreaves. 

[pulseq](https://pulseq.github.io) by University Medical Center Freiburg.


## MATLAB Dependencies
To make sure the code is able to run, the following toolboxes need to be installed: 
[MATLAB](https://www.mathworks.com/products/matlab.html)
[Optimization Toolbox](https://www.mathworks.com/products/optimization.html)
[Signal Processing Toolbox](https://www.mathworks.com/products/signal.html)
[Image Processing Toolbox](https://www.mathworks.com/products/image-processing.html)
[MATLAB Coder](https://www.mathworks.com/products/matlab-coder.html)


 ## Citing
 Zhao Z, Lee NG, Nayak KS. Multidimensional RF Pulse Design with Consideration of Concomitant Field Effects. Magn Reson Med. xxx
 
 ## License
 The Multi_RF packaged code is licensed under a 3-Clause BSD License.

 ## Contact
 If you have any questions, please contact ziweiz@usc.edu

 Ziwei Zhao, University of Southern California, MREL (Magnetic Resonance Engineering Laboratory, PI: Krishna S. Nayak, https://mrel.usc.edu/) May 2024.



