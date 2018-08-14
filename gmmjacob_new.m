function [G,diff_mat]=gmmjacob_new(coef,x,param)

num_dim_x=size(x,2);

% pdf_diff_m=zeros(size(x,1),num_dim_x);
% pdf_diff_s=zeros(size(x,1),num_dim_x);
% for i=1:num_dim_x
% pdf_diff_m(:,i)=-x(:,i).*normpdf(x(:,i),param(i),param(i+num_dim_x));
% pdf_diff_s(:,i)=-(x(:,i).^2-1).*normpdf(x(:,i),param(i),param(i+num_dim_x));
%
% end
diff_mat=zeros(size(x,1),num_dim_x*2);
pdf_total=mvnpdf(x,param(1:num_dim_x),diag(param(num_dim_x+1:2*num_dim_x).^2));
for i=1:num_dim_x
    index=1:num_dim_x;
    index(i)=[];
    
    m=param(i);s=param(i+num_dim_x);
    pdf_diff_m=-(2^(1/2)*exp(-(m - x(:,i)).^2/(2*s^2)).*(2*m - 2*x(:,i)))./(4*pi^(1/2)*s^2*(s^2)^(1/2));
    
    pdf_diff_s= -(2^(1/2)*exp(-(m - x(:,i)).^2/(2*s^2)).*(- m^2 + 2*m*x(:,i) + s^2 - x(:,i).^2))./(2*pi^(1/2)*s*(s^2)^(3/2));
    pdf_cond=pdf_total./normpdf(x(:,i),param(i),param(i+num_dim_x));
%     pdf_cond=mvnpdf(x(:,index),param(index),diag(param(index+num_dim_x).^2));
    diff_mat(:,i)=pdf_cond.*pdf_diff_m;
    diff_mat(:,i+num_dim_x)=pdf_cond.*pdf_diff_s;
end


G=coef'*diff_mat;

end