nataf.coord=zeros(num_smp,num_dim_x);
nataf.mpdf=zeros(num_smp,num_dim_x);

for iteri=1:num_dim_x
    
%   compute marginal pdf of x
%     for itersmp=1:num_smp
%         
%         pdfx(itersmp,iteri)=epipar{itermain,iteri}(ceil(smp_x(itersmp,iteri)/((mesh.mend-mesh.m0)/mesh.m)),2)*smp_x(itersmp,iteri)...
%             +epipar{itermain,iteri}(ceil(smp_x(itersmp,iteri)/((mesh.mend-mesh.m0)/mesh.m)),1);
%         
%     end
%     
    
    pdfx(:,iteri)=epimarpdf(epipar{itermain,iteri},mesh,smp_x(:,iteri));
    
%   compute marginal pdf transformed
    
%     temp=norminv(epicdf(smp_x(:,iteri),epipar{itermain,iteri},mesh));
%     
%     temp(temp==-inf)=-8;temp(isnan(temp))=8;temp(temp==inf)=8;temp(temp<-8)=-8;temp(temp>8)=8;
%     
%     nataf.coord(:,iteri)=temp';
%     
%     clear temp
%     
%     nataf.mpdf(:,iteri)=normpdf(nataf.coord(:,iteri));
end

%---------------------estimatingrho---------------------------------------
nataf.rho=eye(num_dim_x);

% density_rho;

nataf.mvn=mvnpdf(nataf.coord,zeros(1,num_dim_x),nataf.rho);
nataf.mvn(isnan(nataf.mvn))=0;




density_evaluate;