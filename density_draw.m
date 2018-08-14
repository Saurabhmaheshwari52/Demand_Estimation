
% epipar_iteri=a;
function density_draw(mat, mesh, color)
xdim=[];
ydim=[];

for i=1:mesh.m
    xdim = [xdim,mesh.mlist(i):.01:mesh.mlist(i+1)];
    ydim = [ydim,mat(i,2)*(mesh.mlist(i):.01:mesh.mlist(i+1))+mat(i,1)];
end

figure(1)


plot(xdim,ydim,color);
hold on
ylim([0,2.5])
end

    
%         plot(xdim,normpdf(xdim,network.xmean(iteri),sqrt(network.xvar(iteri,iteri))),'k');
%         hold on