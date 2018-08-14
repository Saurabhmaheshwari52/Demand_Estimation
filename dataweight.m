function [weight]=dataweight(y,moment,momentlist)
    n=size(momentlist,1);
    m=size(y,1);
    weight=zeros(n);
    for i=1:n
        for j=1:i
            weight(i,j)=mean((prod(y.^repmat(momentlist(i,:),[m,1]),2)-moment(i)).*(prod(y.^repmat(momentlist(j,:),[m,1]),2)-moment(j)));
            weight(j,i)=weight(i,j);
        end
    end
    
end