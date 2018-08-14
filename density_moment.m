% 1. Data Preparation
num_sample=size(obs_pool_y_old,1)/num_obs;
obs_total_y_old=cell(num_sample,1);
num_dim_y = pca_comps;
for i=1:num_sample     
    obs_total_y_old{i} = obs_pool_y_old((i-1)*num_obs+1:i*num_obs,:);
end
obs_pool_y=zeros(num_obs*num_sample,num_dim_y);
for i=1:num_sample
    obs_total_y{i}=obs_total_y_old{i};
    obs_pool_y((i-1)*num_obs+1:i*num_obs,:)=obs_total_y{i};
end
%-----------------------------------------------------------------------------
% scale sampled y data
smp_y_old_scaled = (smp_y_old - (ones(num_smp,1)*mu_obs))./(ones(num_smp,1)*sigma_obs);
smp_y = smp_y_old_scaled*coeff_pca(:,1:pca_comps);
%-----------------------------------------------------------------------------
% 2. Generate Moments
momentlist_full=zeros(1,num_moment);
for i=1:num_moment
momentlist_temp=nmultichoosek(1:num_dim_y,i);
momentlist_full=[momentlist_full;[momentlist_temp,zeros(size(momentlist_temp,1),num_moment-i)]];
end

mat_temp_index=(1:size(momentlist_full,1))'*ones(1,num_moment);
mat_temp=[momentlist_full(:),mat_temp_index(:)];
momentlist=sparse(mat_temp(mat_temp(:,1)~=0,2),mat_temp(mat_temp(:,1)~=0,1),ones(sum(mat_temp(:,1)~=0),1));

% only first two moments
moment_select=(sum(momentlist>0,2)<=2);
moment_select(1)=0;
%-----------------------------------------------------------------------------
tic
% generate integrated population moment
moment_pdfsmp_y = zeros(size(momentlist,1),1);
pdf_smp = fun_pdf(smp_x(1:num_smp,:));

% Choose numerical sampling method
if ~exist('pdf_aux','var')
    pdf_aux=ones(num_smp,1);
end

weight=1./fun_pdf_smp(smp_x(1:num_smp,:))/num_smp./pdf_aux(1:num_smp,:);
% sample y pdf
for i=1:size(momentlist,1)
    moment_pdfsmp_y(i)=dot(prod(smp_y(1:num_smp,:).^repmat(momentlist(i,:),[num_smp,1]),2).*weight,pdf_smp.*pdf_aux(1:num_smp,:));
end
toc
%-----------------------------------------------------------------------------
% 3. Coefficient Computation
% integrate gx using monte carlo integration
coefpure=zeros(num_smp,size(momentlist,1));
for i=1:size(momentlist,1)

    coefpure(:,i)=prod(smp_y(1:num_smp,:).^repmat(momentlist(i,:),[num_smp,1]),2).*weight;

end 

coefpure=coefpure.*(pdf_aux(1:num_smp,:)*ones(1,size(momentlist,1)));
%-----------------------------------------------------------------------------
tic
% generate observation moment
moment_obs_y = zeros(size(momentlist,1),num_sample);
for iteri=1:num_sample
    obs_y=obs_total_y{iteri};
for i=1:size(momentlist,1)
    moment_obs_y(i,iteri)=mean(prod(obs_y.^repmat(momentlist(i,:),[num_obs,1]),2));
end
end
toc
% generate pool observation moment
moment_obs_pool_y = mean(moment_obs_y,2);
%------------------------------------------------------------------------------
tic
 [sigmadata]=dataweight(obs_pool_y(1:1e4,:),moment_obs_pool_y(moment_select),momentlist(moment_select,:));
toc

sigmasimulation=simulationweight(smp_y(1:1e4,:),1./fun_pdf_smp(smp_x(1:1e4,:)).*fun_pdf(smp_x(1:1e4,:)),moment_obs_pool_y(moment_select),momentlist(moment_select,:));




