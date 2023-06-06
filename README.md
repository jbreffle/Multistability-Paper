# Multistability-Paper
Codes for our paper on "Multistability in neural systems with random cross-connections".

The repository includes separate codes for simulations of finite systems, for mean-field analysis of infinite systems, and for solving approximate analyses of finite systems with binary neural response functions.


Instructions to create paper figures: 

## Figure 1. Single-unit response functions, f(x), produce bistability at s=1
- ./FI_curves.m produces panels A-C

## Figure 2: All activity regimes can be observed in a single random network
- ./example_net_dynamics.m produces panels A-D with the setting "figSim = 'fig2';"
- ./perturb_nets_run_sim.m produces the data for panel E
- ./perturb_nets_plot_sim.m plots the data for panel E

## Figure 3: Smaller networks are multistable at larger g values
- ./simulation-code/sim_rand_nets.m produces the data for panels B-F with appropriate parameter selection
    - Set options as instructed in the "%% Sim options" section
- ./example_net_dynamics.m produces panel A with the setting "figSim = 'fig3';"
- ./same_nets_across_g.m plots the data for panels B-F

## Figure 4: Multistability across phase space
- ./simulation-code/sim_rand_nets.m produces the data for all tanh and logistic unit panels with appropriate parameter selection
    - Set options as instructed in the "%% Sim options" section
- ./simulation-code/sim_rand_nets_binary.m produces the data for all binary unit panels with appropriate parameter selection
        - Set options as instructed in the "%% Sim options" section
- ./grid_sim_plotting.m plots the data produced by the simulation codes

## Supplemental Figure 1: The example network in Fig. 2 is a typical random network
- ./example_net_dynamics.m produces all panels with the setting "figSim = 'fig2';"

## Supplemental Figure 2: The dynamic regimes of individual networks changes as g is scaled
- ./simulation-code/sim_rand_nets.m produces the data with appropriate parameter selection
    - Set options as instructed in the "%% Sim options" section
- ./same_nets_across_g_classes.m plots the data for panels B-F


Notes: 
- The simulation code can run either locally or on a slurm computing cluster with an array of 100 jobs.
- The ./simulation-code/*.sh files run the code on a slurm cluster with each job running one network at each parameter point.
- Data produced by the cluster needs to be combined into a single file with ./process_HPCC_sim.m to be plotted by the analysis scripts
- ./functions/myPlotSettings.m sets the figure settings for each figure
- ./functions/AppendFileNum.m ensures existing data is not overwritten
