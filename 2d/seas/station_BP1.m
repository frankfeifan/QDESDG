clear

addpath ../

data_base_name = '~/scratch/SEAS/BP1_N_4_R_0_P_beta_10_mesh_version_9_dfd436e'; % kozdon.1
% data_base_name = '~/scratch/SEAS/BP1_N_4_R_0_P_beta_10_mesh_version_13_2794eee'; % kozdon.2
% data_base_name = '~/scratch/SEAS/BP1_N_4_R_0_P_beta_10_mesh_version_17_f88d9a3'; % kozdon.3
% data_base_name = '~/scratch/SEAS/BP1_N_4_R_0_P_beta_10_mesh_version_16_0655ba0'; % kozdon.4

load([data_base_name,'_data.mat']);

station_data(yfault, plot_fault{1}, [data_base_name, '_'], ...
            [0, -2.5, -5, -7.5, -10, -12.5, -15, -17.5, -20, -25, -30, -35])
