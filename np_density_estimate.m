%----------------------update each variable------------------------------

for iteri = 1%:num_dim_x
    epiref=epipar{itermain,iteri};
    pdfcond=prod(pdfx(:,[1:iteri-1,iteri+1:end]),2);     
    pdfcond(isnan(pdfcond))=0;
    
    coef=zeros(num_smp,length(momentlist_np(:,1)));
    
    for i=1:length(momentlist_np(:,1))
        coef(:,i)=coefpure_np(:,i).*pdfcond;
    end    
    [epipar_iteri,val_obj(itermain,iteri)] = np_density_optimize(smp_x_subs(:,iteri),...
        moment_y, mesh, coef, weight, pdfcond, epipar{itermain,iteri}, momentlist_np);

    
   %----------------------------------------
    
    pdf_temp = pdfx;
    
    pdf_temp(:,iteri) = epimarpdf(epipar_iteri,mesh,smp_x_subs(:,iteri));
    
    pdf_temp=prod(pdf_temp,2);
    
    moment_pdftemp_y = zeros(size(momentlist_np,1),1);

    for i=1:size(momentlist_np,1)
        moment_pdftemp_y(i)=dot(prod(smp_y_np.^repmat(momentlist_np(i,:),[num_smp,1]),2),pdf_temp)*2^num_dim_x/num_smp;
    end
    
    val_obj_verify(itermain,iteri)= sum(weight.*(moment_pdftemp_y-moment_y).^2)*100;
    
  %-------------------------------------------  
%     epipar_iteri=epipar_iteri./stats(2);
    
   % density_draw;
    
    epipar{itermain+1,iteri}=epipar_iteri;
   
    
%     mse(itermain,iteri)= mse_1d(epipar{itermain,iteri},epipar{itermain+1,iteri},mesh);
    
%     mse_base;
    
%   identifiability test
%     density_mat
    
% % % % % %     density_update
      
end

% density_rho
    