num_smp_k = 1e5;
estimated=1:num_dim_x_old;
num_dim_x=sum(estimated>0); % x dimension
num_param=num_dim_x*2; % number of parameters
if sum(estimated==0)==0
    coef_cpl=coefpure;
else
    coef_cpl = coefpure.*(mvnpdf(smp_x(:,estimated==0),par_mean_x(estimated==0),diag(par_std_x(estimated==0).^2))*ones(1,length(moment_select)));
end

[gmmG_true,diffmat]=gmmjacob_new(coef_cpl(:,moment_select),smp_x(:,estimated>0),[par_mean_x(estimated>0),par_std_x(estimated>0)]);
cov_true=inv(gmmG_true'*((sigmadata+0*num_obs/num_smp_k*sigmasimulation)\gmmG_true))/num_obs;

est_total=zeros(num_sample_t,num_dim_x*2);

tic
options = optimoptions('fmincon','TolFun',1e-12,'MaxFunEvals',5e4,'MaxIter',1e4,'GradObj','on','Display','iter');
est_init=[ones(1,num_dim_x)*1.0,ones(1,num_dim_x)*.2];
fitobj = []; %objective function value
for i=1:2%num_sample_t
    rng(i)
    rand_sim= randsample(num_smp,num_smp_k);
    gmmweight=inv(sigmadata+num_obs/length(rand_sim)*sigmasimulation);
    fitness=@(param)gmmobj(param,smp_x(rand_sim,estimated>0),coef_cpl(rand_sim,moment_select)*num_smp/num_smp_k,moment_obs_y(moment_select,i),gmmweight);
    tic
    [est_total_temp]=  fmincon(fitness,est_init,[],[],[],[] ,[ones(1,num_dim_x)*.5,ones(1,num_dim_x)*0.1],[ones(1,num_dim_x)*1.5,ones(1,num_dim_x)*0.3],[],options);
    est_total(i,:)=est_total_temp;
    fitobj = [fitobj,fitness(est_total(i,:))];
    toc
end
toc


