clear

addpath ../

data_base_name = '~/scratch/SEAS/BP1_N_4_R_0_P_beta_10_mesh_version_5_0e738e1_unclean';
data_base_name = '~/scratch/SEAS/BP1_N_4_R_0_P_beta_10_mesh_version_6_465b541';
data_base_name = '~/scratch/SEAS/BP1_N_4_R_0_P_beta_10_mesh_version_7_465b541';
data_base_name = '~/scratch/SEAS/BP1_N_4_R_0_P_beta_10_mesh_version_8_3ed765a';

data_base_name = '~/scratch/SEAS/BP1_N_4_R_0_P_beta_10_mesh_version_9_dfd436e';
% data_base_name = '~/scratch/SEAS/BP1_N_4_R_0_P_beta_10_mesh_version_10_dfd436e';
% data_base_name = '~/scratch/SEAS/BP1_N_4_R_0_P_beta_10_mesh_version_11_dfd436e';
% data_base_name = '~/scratch/SEAS/BP1_N_4_R_0_P_beta_10_mesh_version_12_dfd436e';

data_base_name = '~/scratch/SEAS/BP1_N_4_R_0_P_beta_10_mesh_version_9_dfd436e';
data_base_name = '~/scratch/SEAS/BP1_N_4_R_1_P_beta_10_mesh_version_9_dfd436e';
data_base_name = '~/scratch/SEAS/BP1_N_8_R_0_P_beta_10_mesh_version_9_dfd436e';




load([data_base_name,'_data.mat']);

hold on
plot_slip_DG(yfault, plot_fault{1}(:), [data_base_name, '_'])
