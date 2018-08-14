function [p]=epimarpdf(epipar,mesh,x)

% p=zeros(size(x));
% 
% for itersmp=1:length(x)
%     
%     p(itersmp)=epipar(ceil(x(itersmp)/((mesh.mend-mesh.m0)/mesh.m)),2)*x(itersmp)...
%         +epipar(ceil(x(itersmp)/((mesh.mend-mesh.m0)/mesh.m)),1);
%     
% end
p = zeros(size(x,1),1);
for i = 1:size(x,1)
    if x(i,1) > mesh.mend || x(i,1) < mesh.m0
        p(i) = 0;
    else
        p(i) = epipar(ceil(x(i,1)/((mesh.mend-mesh.m0)/mesh.m)),2).*x(i,1)+epipar(ceil(x(i,1)/((mesh.mend-mesh.m0)/mesh.m)),1);
    end
end
p(p<0)=0;



