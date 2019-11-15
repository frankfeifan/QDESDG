%{
clear

addpath ../

data_base_name{1} = '~/scratch/SEAS/BP1_N_4_R_0_P_beta_10_mesh_version_9_dfd436e'; % kozdon.1
data_base_name{2} = '~/scratch/SEAS/BP1_N_4_R_0_P_beta_10_mesh_version_13_2794eee'; % kozdon.2
data_base_name{3} = '~/scratch/SEAS/BP1_N_4_R_0_P_beta_10_mesh_version_17_f88d9a3'; % kozdon.3
data_base_name{4} = '~/scratch/SEAS/BP1_N_4_R_0_P_beta_10_mesh_version_16_0655ba0'; % kozdon.4

legend_entry{1} = '800 X 400 :: Vp/2';
legend_entry{2} = '160 X  80 :: Vp/2';
legend_entry{3} = '800 X 400 :: free';
legend_entry{4} = '160 X  80 :: free';

figure(1)
for k = 1:length(data_base_name)
  disp(legend_entry{k})
  compute_event_time([data_base_name{k}, '_'], [180, 230])
  hold on
end
hold off

legend(legend_entry)

ylabel('max(|V|')
xlabel('time (years)')
axis tight
xlim([180 230])
% matlab2tikz('free_Vp_compare.tikz');
%}

clear

data_base_name{1} = '~/scratch/SEAS/BP1_N_4_R_0_P_beta_10_mesh_version_9_dfd436e'; % kozdon.1
data_base_name{2} = '~/scratch/SEAS/BP1_N_4_R_0_P_beta_10_mesh_version_13_2794eee'; % kozdon.2
data_base_name{3} = '~/scratch/SEAS/BP1_N_4_R_0_P_beta_10_mesh_version_18_eaa88a4';
data_base_name{4} = '~/scratch/SEAS/BP1_N_4_R_0_P_beta_10_mesh_version_19_eaa88a4';

legend_entry{1} = '800 X 400 :: Vp/2';
legend_entry{2} = '160 X  80 :: Vp/2';
legend_entry{3} = '800 X  80 :: Vp/2';
legend_entry{4} = '160 X 400 :: Vp/2';

figure(2)
for k = 1:length(data_base_name)
  disp(legend_entry{k})
  compute_event_time([data_base_name{k}, '_'], [180, 230])
  hold on
end
hold off

legend(legend_entry)

ylabel('max(|V|')
xlabel('time (years)')
axis tight
xlim([180 230])
% matlab2tikz('Vp_compare.tikz');

clear

data_base_name{3} = '~/scratch/SEAS/BP1_N_4_R_0_P_beta_10_mesh_version_17_f88d9a3'; % kozdon.3
data_base_name{4} = '~/scratch/SEAS/BP1_N_4_R_0_P_beta_10_mesh_version_16_0655ba0'; % kozdon.4
data_base_name{1} = '~/scratch/SEAS/BP1_N_4_R_0_P_beta_10_mesh_version_20_eaa88a4';
data_base_name{2} = '~/scratch/SEAS/BP1_N_4_R_0_P_beta_10_mesh_version_21_eaa88a4';

legend_entry{3} = '800 X 400 :: free';
legend_entry{4} = '160 X  80 :: free';
legend_entry{1} = '800 X  80 :: free';
legend_entry{2} = '160 X 400 :: free';

figure(3)
for k = 1:length(data_base_name)
  disp(legend_entry{k})
  compute_event_time([data_base_name{k}, '_'], [180, 230])
  hold on
end
hold off

legend(legend_entry)

ylabel('max(|V|')
xlabel('time (years)')
axis tight
xlim([180 230])
% matlab2tikz('free_compare.tikz');
