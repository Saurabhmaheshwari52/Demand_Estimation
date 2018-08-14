
pdf_est=ones(num_smp,1);

for iteri=1:num_dim_x
    
    pdf_est=pdf_est.*pdfx(:,iteri)./nataf.mpdf(:,iteri);
    
end


pdf_est=pdf_est.*nataf.mvn;

pdf_est(isnan(pdf_est))=0;

val_msetot(itermain)=mean((pdf_est-pdf_smp).^2.*pdf_smp);

val_mseshf(itermain)=mean((pdf_est-pdf_old).^2.*pdf_smp);
    
val_msemar(itermain,:)=mean((pdfx-pdf_mar).^2.*pdf_mar);


pdf_old=pdf_est;
