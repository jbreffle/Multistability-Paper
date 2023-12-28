# Multistability-Paper

## Description

This repository contains all code necessary to produce the figures of the paper:

```text
Breffle J, Mokashe S, Qiu S, Miller P. 2023. Multistability in neural systems with random cross-connections. Biol Cybern. 117(6):485–506. doi:10.1007/s00422-023-00981-w.
```

Journal link: <https://doi.org/10.1007/s00422-023-00981-w> \
Link to publicly available view-only version: <https://rdcu.be/dujJP>

The repository includes separate codes for simulations of finite systems, for mean-field analysis of infinite systems, and for solving approximate analyses of finite systems with binary neural response functions.

## Figure-by-figure instructions

### Figure 1: Single-unit response functions, $f(x)$, produce bistability at $s=1$

- ```./figure-code/FI_curves.m``` produces panels A-C

### Figure 2: All activity regimes can be observed in a single random network

- ```./figure-code/example_net_dynamics.m``` produces panels A-D with the setting "figSim = 'fig2';"
- ```./figure-code/perturb_nets_run_sim.m``` produces the data for panel E
- ```./figure-code/perturb_nets_plot_sim.m``` plots the data for panel E

### Figure 3: Smaller networks are multistable at larger $g$ values

- ```./figure-code/simulation-code/sim_rand_nets.m``` produces the data for panels B-F with appropriate parameter selection
  - Set options as instructed in the "%% Sim options" section
- ```./figure-code/example_net_dynamics.m``` produces panel A with the setting "figSim = 'fig3';"
- ```./figure-code/same_nets_across_g.m``` plots the data for panels B-F

### Figure 4: Multistability across phase space

- ```./figure-code/simulation-code/sim_rand_nets.m``` produces the data for all tanh and logistic unit panels with appropriate parameter selection
  - Set options as instructed in the "%% Sim options" section
- ```./figure-code/simulation-code/sim_rand_nets_binary.m``` produces the data for all binary unit panels with appropriate parameter selection
        - Set options as instructed in the "%% Sim options" section
- ```./figure-code/grid_sim_plotting.m``` plots the data produced by the simulation codes

### Figure 5: Zipf's Law from data and from sampling a log-normal distribution

- all codes are found in the ```./Zipf_Figure_Codes/``` folder

### Figure 6: Infinite limit calculations with logistic f-I curves

- all codes are found in the ```./Inf_Logistic_Codes/``` folder
- the code ```inf_logistic_grid.m``` should run if the other functions are in the same directory
- the code ```plot_phase_diagram.m``` should produce the figure so long as its input file is adjusted to match the output file of ```inf_logistic_grid.m```

### Figure 8: Probability of multistability in finite circuits

- ```erf_code_binary_units1.m``` produces continuous analytic curves
- ```binary_net_sim_count.m``` produces simulated data points

### Figure 9: Infinite limit calculations for binary response curves

- ```erfApprox_maxk_infN.m``` produces the analytic results
- the binary simulation code of Figure 5 or 8 can be adapted to produce Figure 9C

### Supplemental Figure 1: The example network in Fig. 2 is a typical random network

- ```./figure-code/example_net_dynamics.m``` produces all panels with the setting ```figSim='fig2';```

### Supplemental Figure 2: The dynamic regimes of individual networks changes as g is scaled

- ```./figure-code/simulation-code/sim_rand_nets.m``` produces the data with appropriate parameter selection
  - Set options as instructed in the "%% Sim options" section
- ```./figure-code/same_nets_across_g_classes.m``` plots the data for panels B-F

## Usage notes

- The simulation code can run either locally or on a slurm computing cluster with an array of 100 jobs.
- The ```./figure-code/simulation-code/*.sh``` files run the code on a slurm cluster with each job running one network at each parameter point.
- Data produced by the cluster needs to be combined into a single file with ```./figure-code/process_HPCC_sim.m``` to be plotted by the analysis scripts
- ```./figure-code/functions/myPlotSettings.m``` sets the figure settings for each figure
- ```./figure-code/functions/AppendFileNum.m``` ensures existing data is not overwritten
