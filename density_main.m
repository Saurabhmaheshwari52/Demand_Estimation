% multivariate probability density estimation based on moment matching
% main function

% define mapping function
num_dim_x = 3;
num_dim_y = 3;
% fun_map=@(x) x.^2;
%fun_map = @(x) [x.^2/2,x(:,1).^3/4+x(:,2).^3/4,x(:,3).^3/4+x(:,4).^3/4];
fun_map = @(x) braess(x,num_dim_y);

% define true pdf
density_true;


% draw sample
rng(1);   
cur_seed=rng;
num_smp = 1e3;
smp_x=[];smp_y=[];
for i=1:1
    rng(cur_seed);
    smp_x_temp = unifrnd(0,2,num_smp,num_dim_x);
    smp_y_temp = fun_map(smp_x_temp);
    smp_x=[smp_x;smp_x_temp];
    smp_y=[smp_y;smp_y_temp];
    cur_seed=rng;
end
% save eval_braess

% load eval_braess



% define reference marginal pdf
fun_ref = @(x)ones(size(x))/(2);

% synthesize data
rng(0);    
num_obs = 1e3;
obs_x = mvnrnd(par_mean_x,par_var_x,num_obs);
obs_y = fun_map(obs_x);

% generate moment_list
global momentlist
momentlist=[];
num_moment = 2;
for i=0:num_moment
    enummoment(1,[],i,num_dim_y);
end


%%%%%% Aren't moment_pdfobs_y and moment_obs_y the same thing?

% generate observed y moment with known pdf
moment_pdfobs_y = zeros(size(momentlist,1),1);
% pdf_obs = fun_pdf(obs_x,par_mean_x,par_var_x);
pdf_obs= ones(num_obs,1)*1/num_obs;
for i=1:size(momentlist,1)
    moment_pdfobs_y(i)=dot(prod(obs_y.^repmat(momentlist(i,:),[num_obs,1]),2),pdf_obs);
end 

% generate observation moment
moment_obs_y = zeros(size(momentlist,1),1);
for i=1:size(momentlist,1)
    moment_obs_y(i)=mean(prod(obs_y.^repmat(momentlist(i,:),[num_obs,1]),2));
end 

%%% parameters in density_mesh?
% setting epi spline with Gaussian quadrature points
density_mesh;



% initialize epispline
fun_mar=@(x)normpdf(x,par_mean_x(1),sqrt(par_var_x(1,1))); % marginal density
% 
a=density_epiapprox(fun_mar,mesh); % the best we can do
%%% what is a?

% generate sampled y moment with known pdf
num_smp=length(smp_y(:,1));
moment_pdfsmp_y = zeros(size(momentlist,1),1);
pdf_smp = fun_pdf(smp_x);
for i=1:size(momentlist,1)
    moment_pdfsmp_y(i)=dot(prod(smp_y.^repmat(momentlist(i,:),[num_smp,1]),2),...
    pdf_smp)*2^num_dim_x/num_smp;
end

%%%%%%% why 2^num_dim_x/num_smp??

pdf_mar=zeros(size(smp_x));
for i=1:num_dim_x
    pdf_mar(:,i) = normpdf(smp_x(:,i),par_mean_x(i),sqrt(par_var_x(i,i)));
end



% integrate gx using monte carlo integration

for i=1:size(momentlist,1)

    coefpure(:,i)=prod(smp_y.^repmat(momentlist(i,:),[num_smp,1]),2)*2^num_dim_x/num_smp;

end 

% set initial density to be uniform distributed
for iteri = 1:num_dim_x
    epipar{1,iteri}= [ones(mesh.m,1)*(1/(mesh.mend-mesh.m0)),zeros(mesh.m,1)];
end
pdf_old=1/(mesh.mend-mesh.m0)^num_dim_x;


% moment matching weight
weight=1.^sum(momentlist,2);

% decomposition algorithm
moment_y=moment_pdfsmp_y;

num_itermain = 1000;
val_obj = zeros(num_itermain,num_dim_x);
val_msemar = zeros(num_itermain,num_dim_x);

% Gauss Seidel need to turn off density update in estimate

% for itermain = 1:num_itermain
%     density_construct; density_estimate;
% end density_construct;

% Gauss Jordan need to turn on density update in estimate
itermain = 1;
density_construct
for itermain = 1:num_itermain
    itermain
    density_estimate; 
    density_evaluate;
    if itermain>1 && val_obj(itermain-1,1)-val_obj(itermain,1)<1e-5
        break
    end
end

% plotting
density_draw;

