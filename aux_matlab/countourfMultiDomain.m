function countourfMultiDomain(mesh,data,ncontour)
    %plots data on multiDomains meshes
    hold on;
    if ~exist('ncontour'); ncontour=32;end
    
    minC=min(data(:));
    maxC=max(data(:));
    contours=linspace(minC,maxC,ncontour);
    for i=1:mesh.domains.n
        X=mesh.domains.X{i};
        Y=mesh.domains.Y{i};
        C=zeros(size(X));
        C(:)=data(mesh.domains.indexes{i});
        contourf(X,Y,C,contours,'linecolor','none');
        caxis([minC,maxC+1e-9]);
    end
