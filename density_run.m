
load friedrichshain_od.mat

dest=repmat([1:23],1,23);
orig=[];
for i=1:23
    orig=[orig,ones(1,23)*i];
end
index=(orig==dest);
orig=orig(index==0);
dest=dest(index==0);
%%%%%%%% To decrease the number of unknown OD pairs
% % % [sortedOD, I] = sort(mean_x_true, 'descend');
% % % x_threshold = sortedOD(21,1);
% % % mean_x_true80 = mean_x_true(mean_x_true > 40);
% % % [sortedOD80,I80] = sort(mean_x_true80, 'descend');
%%%%%%%%
x_threshold=40;
y_threshold=500;

estimated_old=1:23^2;
estimated_old=estimated_old(index==0);
estimated_old=estimated_old(mean_x_true>x_threshold);

orig=orig(mean_x_true>x_threshold);
dest=dest(mean_x_true>x_threshold);

num_dim_x_old = sum(mean_x_true>x_threshold);
num_dim_y_old = sum(mean_y_true>y_threshold);
observed_old=1:sum(mean_y_true>-1e-16);
observed_old(mean_y_true<(y_threshold+1e-16))=[];

mean_x_true=mean_x_true(mean_x_true>x_threshold)';
mean_y_true=mean_y_true(mean_y_true>y_threshold)';

load friedrichshain_sim_500_40.mat
num_smp=1e5;

smp_y_old=smp_y_old(1:num_smp,:);
smp_x_old=smp_x_old(1:num_smp,:);
% % % smp_x_old = smp_x_old(1:num_smp,I80(1:20,:));
smp_x=smp_x_old./(ones(num_smp,1)*mean_x_true);

% 3. define true pdf
density_model;

% 4. load data sample

num_obs = 250;

num_sample_t=200;

record_i=0;
est_record=cell(9,1);
gmm_record=cell(9,1);
omega_record=cell(9,1);

load friedrichshain_obs_test1.mat
obs_pool_y_old=obs_pool_y_old(1:num_obs*num_sample_t,:);
[zscore_obs, mu_obs, sigma_obs] = zscore(obs_pool_y_old);
[coeff_pca,obs_pool_y_pca, evs] = princomp(zscore_obs);
evs_std = cumsum(evs)/sum(evs);
%which99 = min(find(evs_std > 0.99));
which99 = 17
obs_pool_y_old = obs_pool_y_pca(:,1:which99);

for sensor_i=[750]
    for moment_i= 2
%         observed=1:num_dim_y_old;
%         observed(mean_y_true<sensor_i)=0;
%         [~,idx_temp]=licols(smp_y_old(1:1e4,observed>0),1e-8);
%         observed_temp=observed(observed>0);
%         observed=zeros(1,num_dim_y_old);
%         observed(observed_temp(idx_temp))=observed_temp(idx_temp);
%         observed_old(observed==0)=[];
%         
        num_moment=moment_i;
        density_moment;
        density_estimate;
        record_i=record_i+1;
        
        
        omega_record{record_i}=inv(gmmG_true'*inv(sigmadata+num_obs/num_smp_k*sigmasimulation)*gmmG_true);
    end
end
