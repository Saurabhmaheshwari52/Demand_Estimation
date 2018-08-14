par_mean_x = ones(1,num_dim_x_old);
% mean_x_change=40;
% mean_x_change = 0;
% par_mean_x(mean_x_true>mean_x_change & mod(orig+dest,2)==0)=.9;
% par_mean_x(mean_x_true>mean_x_change & mod(orig+dest,2)==1)=1.1;
% 
% par_mean_x(mean_x_true>mean_x_change)=.9+0.2/(sum(mean_x_true>mean_x_change)-1)*(0:sum(mean_x_true>mean_x_change)-1)';
par_std_x=par_mean_x*.2;
 
% par_mean_x(mean_x_true>mean_x_change)=.9+(orig(mean_x_true>mean_x_change)-1)/22*.2;
% par_std_x(mean_x_true>mean_x_change)=.18+(dest(mean_x_true>mean_x_change)-1)/22*.04;


par_var_x = diag((par_std_x).^2);

fun_pdf= @(x)mvnpdf(x,par_mean_x,par_var_x);

fun_pdf_smp=@(x)mvnpdf(x,ones(1,num_dim_x_old),diag(0.04*ones(num_dim_x_old,1)));