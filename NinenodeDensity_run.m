%clear
load('ninenodes_integrate_naive.mat')
load('ninenodes_flownormal.mat')
load('mesh.mat')
mean_x_true = mean(smp_x_old);
mean_y_true = mean(smp_y_old);
x_threshold=0;
y_threshold=0;
num_dim_x_old = sum(mean_x_true>x_threshold);
num_dim_y_old = sum(mean_y_true>y_threshold);
num_smp=1e6;
smp_y_old=smp_y_old(1:num_smp,:);
smp_x_old = smp_x_old(1:num_smp,:);
smp_x=smp_x_old./(ones(num_smp,1)*mean_x_true); % scaling sample x
%-----------------------------------------------------------------------------
% 3. define true pdf
density_model;
num_obs = 2500;
num_sample_t=200;
record_i=0;
est_record=cell(9,1);
gmm_record=cell(9,1);
omega_record=cell(9,1);
%-----------------------------------------------------------------------------
% PCA 
%y_indices = [1,2,3,4,5,7,9,12,13,14];
y_indices = [1,5,7,9,12,4,13];
obs_pool_y_old=flow_normal_total(1:num_obs*num_sample_t,y_indices);
[zscore_obs, mu_obs, sigma_obs] = zscore(obs_pool_y_old);
[coeff_pca,obs_pool_y_pca, evs] = princomp(zscore_obs);
pca_comps = 6; % first 4 components
obs_pool_y_old = obs_pool_y_pca(:,1:pca_comps);
%-----------------------------------------------------------------------------
% parametric estimation
for moment_i = 4
        num_moment=moment_i;
        density_moment;
        par_density_estimate;
        record_i=record_i+1;
        omega_record{record_i}=inv(gmmG_true'*inv(sigmadata+num_obs/num_smp_k*sigmasimulation)*gmmG_true);
end
% non-parametric estimation
np_density_main;
