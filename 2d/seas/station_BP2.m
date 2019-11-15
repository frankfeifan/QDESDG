clear

addpath ../

files = {...
'~/scratch/qdes/BP2/BP2_N_4_R_0_P_beta_10_mesh_version_1_f8ee70b_unclean',
'~/scratch/qdes/BP2/BP2_N_4_R_1_P_beta_10_mesh_version_1_f8ee70b_unclean',
'~/scratch/qdes/BP2/BP2_N_4_R_2_P_beta_10_mesh_version_1_f8ee70b_unclean',
'~/scratch/qdes/BP2/BP2_N_4_R_3_P_beta_10_mesh_version_1_f8ee70b_unclean',
'~/scratch/qdes/BP2/BP2_N_4_R_4_P_beta_10_mesh_version_1_f8ee70b_unclean',
'~/scratch/qdes/BP2/BP2_N_4_R_0_P_beta_10_mesh_version_2_f8ee70b_unclean',
'~/scratch/qdes/BP2/BP2_N_4_R_0_P_beta_10_mesh_version_3_f8ee70b_unclean',
'~/scratch/qdes/BP2/BP2_N_4_R_0_P_beta_10_mesh_version_4_f8ee70b_unclean',
'~/scratch/qdes/BP2/BP2_N_2_R_1_P_beta_10_mesh_version_2_f8ee70b_unclean',
'~/scratch/qdes/BP2/BP2_N_2_R_2_P_beta_10_mesh_version_2_f8ee70b_unclean',
'~/scratch/qdes/BP2/BP2_N_2_R_3_P_beta_10_mesh_version_2_f8ee70b_unclean',
'~/scratch/qdes/BP2/BP2_N_2_R_4_P_beta_10_mesh_version_2_f8ee70b_unclean',
'~/scratch/qdes/BP2/BP2_N_2_R_0_P_beta_10_mesh_version_3_f8ee70b_unclean',
'~/scratch/qdes/BP2/BP2_N_2_R_0_P_beta_10_mesh_version_4_f8ee70b_unclean',
'~/scratch/qdes/BP2/BP2_N_2_R_0_P_beta_10_mesh_version_5_f8ee70b_unclean',
}

for k = 1:length(files)

  data_base_name = files{k};

  load([data_base_name,'_data.mat']);

  station_data(yfault, plot_fault{1}, [data_base_name, '_'], ...
  [0, -2.4, -4.8, -7.2, -9.6, -12, -14.4, -16.8, -19.2, -24, -28.8, -36])
end
