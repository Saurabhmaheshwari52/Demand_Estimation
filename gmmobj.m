function [val,grad] = gmmobj(param,x,coef,moment,weight)

num_dim_x=size(x,2);
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
diff=(mvnpdf(x,param(1:num_dim_x),diag(param(num_dim_x+1:2*num_dim_x).^2))'*coef)'-moment;
% Braess
% diff=((ones(size(x))*1/param.*(x<=param))'*coef)'-moment;
val=diff'*weight*diff;

grad=2*diff'*weight*gmmjacob_new(coef,x,param);
end

