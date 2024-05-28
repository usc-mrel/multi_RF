# Multidimensional RF Pulse Design with Consideration of Concomitant Field Effects
 
Code for generating simulation results for the original and proposed methods in this paper: *Ziwei Zhao, Nam G. Lee, Krishna S. Nayak. "Multidimensional RF Pulse Design with Consideration of Concomitant Field Effects"*. Submitted to Mag. Reson. Med. 2024

## Set up

### Initialization
 - Please edit the path of the third party packages in **sim_concomitantfield_approx.m** line 11.
 - Please edit the path of the third party packages in **sim_concomitantfields_8channel.m** line 10-14.
 - Please edit the path of the third party packages in **sim_concomitantfields_1channel.m** line 13-19.
 - Please changes the path of the b0 and b1 maps in **sim_concomitantfields_8channel.m** line 48.
 - Please changes the path of the GIRF estimation in **sim_concomitantfields_8channel.m** line 82.
 - Please uncomment the code for selecting mask, b0 and b1 maps at different isocenters in **sim_concomitantfields_1channel.m** line 47-95.
 - Please check line 103 in **sim_concomitantfields_1channel.m** for appropriate off isocenter settings to match with relating mask, b1 and b0 map selections. 

## Code Structure
 
### Figure script: 
- Figure 2 was generated using **./sim/sim_concomitantfield_approx.m**.
- Figure 3 was generated using **./sim/sim_concomitantfields_8channel.m** as one sample case. 
- Figure 4 was generated using **./figures/field_strength/mrm_figure3_fieldstrength.m**. The relating *.mat* files were generated using **./sim/sim_concomitantfields_8channel.m**
- Figure 5 (A) simulation part was generated using **./figures/sim_results/phantom_055T/mrm_figure4_tailored_b0b1.mat**. The relating *.mat* files were generated using **./sim/sim_concomitantfields_1channel.m**
- Supporting Figure S1 was generated using **./sim_results/distance/mrm_sup_figure1_distance.m**
- Supporting Figure S2 was generated using **./sim_results/duration/mrm_sup_figure2_duration_v1.m**

### Demo script: 
- **sim_concomitantfield_approx.m** generates concomitant field accuracy results (Figure 2).
- **sim_concomitantfields_8channel.m** computes the excitation profiles using original and proposed methods at different isocenter and different main field strengths using 8 channel setups. The results can be saved as *.mat* format.
- **sim_concomitantfields_1channel.m** computes the excitation profiles using original and proposed methods at different isocenter and different main field strengths using 1 channel setups. The results can be saved as *.mat* format.

### B0, B1+ and GIRF measurement datasets: 
- GIRF measurements at 0.55T are provided in **./b0b1_map_055T/GIRF_20200221_Duyn_method_coil2.mat** 
- B0, B1+ maps with masks at 0.55T with different isocenter (0, 5, 10, 15cm) are provided in **./b0b1_map_055T** folder.
- **./b0b1_map_055T/gene_b0_map_flash3d.mat** is used to estimate a B0 map from flash_2d images using linear fitting method.

### 2D Spin echo sequence:
- **./seq/demo_pulseq_Pauly_1989_JMR_modified.mat** generates the 2D single slice spin echo sequence in Pulseq *.seq* file with the excitation RF replaced with the designed RF pulse.

### Functions: 
- **./sim/funcs/STA_maxwell_system_matrix_con.m** calculates a small-tip-angle system matrix with the consideration of concomitant fields during the iterative RF design (proposed method).
- **./sim/funcs/STA_maxwell_system_matrix.m** calculates a small-tip-angle system matrix without the consideration of concomitant fields during the iterative RF design (original method).
- **./sim/funcs/figure_out_transformation_matrix.m** and **./sim/funcs/geometry** folder outputs the coordination transformation matrix that is used in the application of GIRF at 0.55T.
- **./sim/funcs/apply_GIRF_tx.m** calculates the GIRF-corrected trajectory.
- **./sim/funcs/calcgradinfo.m** calculates the kspace trajectory from the input gradient waveforms.
- **./seq/input_pulseq_cartesian_spin_echo_multi_RF.m** saves sequence parameters that is used in the pulseq file.
- **./seq/calculate_cartesian_spin_echo_imaging_parameters.m** calculates basic sequence parameters and define gradient performance at 0.55T that are used in the pulseq file.


## Dependencies
[reVERSE_GIRF](https://github.com/mriphysics/reverse-GIRF?tab=readme-ov-file) by Shaihan Malik, July 2015.

[Phase_relaxed_CPMG_excitation](https://github.com/mriphysics/phase_relaxed_CPMG_excitation) by Shaihan Malik, July 2015.

[lsqrSOL](https://github.com/areslp/matlab/tree/master/lsqrSOL) for solving a linear non-square matrix problem. 

[Bloch_simulation](http://mrsrl.stanford.edu/~brian/blochsim/) by Brian Hargreaves. 


 ## Citing
 Zhao Z, Lee NG, Nayak KS. Multidimensional RF Pulse Design with Consideration of Concomitant Field Effects. Magn Reson Med. xxx
 
 ## License
 The Improved_3DRT_Speech packaged code is licensed under a 3-Clause BSD License.

 ## Contact
 If you have any questions, please contact ziweiz@usc.edu

 Ziwei Zhao, University of Southern California, MREL (Magnetic Resonance Engineering Laboratory, PI: Krishna S. Nayak, https://mrel.usc.edu/) May 2024.



