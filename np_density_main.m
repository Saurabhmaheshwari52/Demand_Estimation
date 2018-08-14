% multivariate probability density estimation based on moment matching
% main function

num_dim_x = 4;
%------generate y obs data-----------------------------------------------------------
num_obs = 2500;
num_sample_t=size(flow_normal_total,1)/num_obs;
num_dim_y_np = 10; % 6 components for np estimation
obs_pool_y_old_np=flow_normal_total(1:num_obs*num_sample_t,y_indices);
obs_total_y_old_np=cell(num_sample_t,1);
for i=1:num_sample_t
obs_total_y_old_np{i} = obs_pool_y_old_np((i-1)*num_obs+1:i*num_obs,:);
end
for i=1:num_sample_t
obs_total_y_np{i}=obs_total_y_old_np{i}./(ones(num_obs,1)*mean_y_true(y_indices))-1;
obs_pool_y_np((i-1)*num_obs+1:i*num_obs,:)=obs_total_y_np{i};
end
% obs_pool_y_old_np = obs_pool_y_pca(:,1:num_dim_y_np);
% obs_y_np = obs_pool_y;
%-------------------------------------------------------------------------------------
% generate moments
num_moment = 5;
momentlist_full=zeros(1,num_moment);
for i=1:num_moment
momentlist_temp=nmultichoosek(1:num_dim_y_np,i);
momentlist_full=[momentlist_full;[momentlist_temp,zeros(size(momentlist_temp,1),num_moment-i)]];
end

mat_temp_index=(1:size(momentlist_full,1))'*ones(1,num_moment);
mat_temp=[momentlist_full(:),mat_temp_index(:)];
momentlist_np=sparse(mat_temp(mat_temp(:,1)~=0,2),mat_temp(mat_temp(:,1)~=0,1),ones(sum(mat_temp(:,1)~=0),1));

% only first four moments
moment_select=(sum(momentlist_np>0,2)<=num_moment);
moment_select(1)=0;
%---------------------------------------------------------------------------------
% generate observed y moment with known pdf
% load('observed_moment_8.mat')
moment_obs_y_np = zeros(size(momentlist_np,1),num_sample_t);
for iteri=1:num_sample_t
    obs_y_np=obs_total_y_np{iteri};
for i=1:size(momentlist_np,1)
    moment_obs_y_np(i,iteri)=mean(prod(obs_y_np.^repmat(momentlist_np(i,:),[num_obs,1]),2));
end
end
%---------------------------------------------------------------------------------
% density_mesh;
load('mesh.mat');
%---------------------------------------------------------------------------------
% generate sampled y moment with known pdf
num = 10000;
rng(1)
samples = randsample(size(smp_y_old,1), num);
smp_y_old_subs = smp_y_old(samples,y_indices);
smp_x_subs = smp_x(samples,:);
num_smp=length(smp_y_old_subs(:,1));
smp_y_np=smp_y_old_subs(1:num_smp,:)./(ones(num_smp,1)*mean_y_true(y_indices))-1;
% smp_y_old_scaled_np = (smp_y_old_subs - (ones(num_smp,1)*mu_obs))./(ones(num_smp,1)*sigma_obs);
% smp_y_np = smp_y_old_scaled_np*coeff_pca(:,1:num_dim_y_np);
moment_pdfsmp_y_np = zeros(size(momentlist_np,1),1);
pdf_smp = fun_pdf(smp_x_subs(1:num_smp,:));
pdf_aux=ones(num_smp,1);
weight_pdf=1./fun_pdf_smp(smp_x_subs(1:num_smp,:))/num_smp./pdf_aux(1:num_smp,:);
for i=1:size(momentlist_np,1)
    moment_pdfsmp_y_np(i)=dot(prod(smp_y_np(1:num_smp,:).^repmat(momentlist_np(i,:),[num_smp,1]),2).*weight_pdf,pdf_smp.*pdf_aux(1:num_smp,:));
end
%---------------------------------------------------------------------------------
% integrate gx using monte carlo integration
coefpure_np=zeros(num_smp,size(momentlist_np,1));
for i=1:size(momentlist_np,1)
    coefpure_np(:,i)=prod(smp_y_np(1:num_smp,:).^repmat(momentlist_np(i,:),[num_smp,1]),2).*weight_pdf;
end 
coefpure_np=coefpure_np.*(pdf_aux(1:num_smp,:)*ones(1,size(momentlist_np,1)));
%---------------------------------------------------------------------------------
% moment matching weight
weight=1.^sum(momentlist_np,2);

fun_mar=@(x,mean,std)normpdf(x,mean,std); % marginal density
pdfx = zeros(num_smp, num_dim_x);
for i = 1:num_dim_x
        pdfx(:,i) = fun_mar(smp_x_subs(:,i),par_mean_x(i),par_std_x(i));
end
%---------------------------------------------------------------------------------
% set initial density to be uniform distributed
epipar = {};
for iteri = 1:num_dim_x
    epipar{1,iteri}= [ones(mesh.m,1)*(1/(mesh.mend-mesh.m0)),zeros(mesh.m,1)];
end
pdf_old=1/(mesh.mend-mesh.m0)^num_dim_x;
coeffs = cell(num_sample_t,1);
for sp = 1:5%num_sample_t
    sp
    % decomposition algorithm
%     moment_y=moment_pdfsmp_y_np(:,sp);
%     pd = makedist('Normal');
%     t = truncate(pd,0,2);
%     moment_y = random(t,size(moment_obs_y_np,1),1);
% binornd(1, 0.5, 3003,1);
    moment_y = moment_pdfsmp_y_np(:,sp);%.*random(t,size(moment_obs_y_np,1),1);
%     maxi = find(moment_y == max(moment_y));
%     moment_y(1717) = 10^8;
    num_itermain = 1000;
    val_obj = zeros(num_itermain,num_dim_x);
    val_msemar = zeros(num_itermain,num_dim_x);
    for itermain = 1:num_itermain
        np_density_estimate; 
%       density_evaluate;
        if itermain>1 && val_obj(itermain-1,1)-val_obj(itermain,1)<1e-5
            break
        end
    end
    coeffs{sp,1} = epipar_iteri;
end
% plotting density for the last sample
density_draw(coeffs{2,1}, mesh, 'green')
mat = [smp_x_subs(:,1),pdfx(:,1)];
mat = sortrows(mat, 1);
plot(mat(:,1),mat(:,2),'-r')
%------------------------------------------------------------------------------
% generating obs moments for epi spline
inter = mesh.mlist;
epi_moments = zeros(4,num_sample_t);
for j = 1:4
    differ = [inter.^(j+1)/(j+1); inter.^(j+2)/(j+2)];
    differ = diff(differ,1,2);
    for i = 1:num_sample_t
        epi_moments(j,i)=sum(sum(differ.*coeffs{i,1}'));
    end
end

%------------------------------------------------------------------------------
figure(1)
plot(1:3003, moment_pdfsmp_y_np(:,1), 'or')
hold on
figure(1)
plot(1:3003, moment_obs_y_np(:,1), 'ob')