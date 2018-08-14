%----------------------update each variable------------------------------

for iteri = 1:num_dim_x
    
    
    
    epiref=epipar{itermain,iteri};
    
 % independent case  
    pdfcond=prod(pdfx(:,[1:iteri-1,iteri+1:end]),2);
 % dependent case
%     pdfcond=prod(pdfx(:,[1:iteri-1,iteri+1:end]),2)./prod(nataf.mpdf(:,[1:iteri-1,iteri+1:end]),2).*nataf.mvn;
% dependent case
%     pdfcond=ones(num_smp,1);
%     
%     
%     
%     for iterj=1:num_dim_x
%         
%         if iterj==iteri
%             
%         else
%             
%             pdfcond=pdfcond.*pdfx(:,iterj)./nataf.mpdf(:,iterj);
%             
%         end
%         
%     end
%     
%     pdfcond=pdfcond./nataf.mpdf(:,iteri).*nataf.mvn;
    
% independence case
%     for iterj=1:num_dim_x
%         
%         if iterj==iteri
%             
%         else
%             
%             pdfcond=pdfcond.*pdfx(:,iterj);
%             
%         end
%         
%     end
      
    pdfcond(isnan(pdfcond))=0;
    
    coef=zeros(num_smp,length(momentlist(:,1)));
    
    for i=1:length(momentlist(:,1))
        coef(:,i)=coefpure_np(:,i).*pdfcond;
    end
    
%     coeffi{itermain,iteri}=coef;

 % compute coefficient  
    
%     [epipar_iteri0, stats0, moment_est0] = epispline_poly(smp_x(:,iteri), moment_y, mesh, coef, weight, epiref);
    
    [epipar_iteri,val_obj(itermain,iteri)] = density_optimize(smp_x(:,iteri), moment_y, mesh, coef, weight, pdfcond, epipar{itermain,iteri});

    
   %----------------------------------------
    
    pdf_temp = pdfx;
    
    pdf_temp(:,iteri) = epimarpdf(epipar_iteri,mesh,smp_x(:,iteri));
    
    pdf_temp=prod(pdf_temp,2);
    
    moment_pdftemp_y = zeros(size(momentlist,1),1);

    for i=1:size(momentlist,1)
        moment_pdftemp_y(i)=dot(prod(smp_y.^repmat(momentlist(i,:),[num_smp,1]),2),pdf_temp)*2^num_dim_x/num_smp;
    end
    
    val_obj_verify(itermain,iteri)= sum(weight.*(moment_pdftemp_y-moment_y).^2)*100;
    
  %-------------------------------------------  
%     epipar_iteri=epipar_iteri./stats(2);
    
%     density_draw;
    
    epipar{itermain+1,iteri}=epipar_iteri;
    
    mse(itermain,iteri)= mse_1d(epipar{itermain,iteri},epipar{itermain+1,iteri},mesh);
    
%     mse_base;
    
%   identifiability test
%     density_mat
    
    density_update
      
end

% density_rho
    