function mesh = CreateLMesh(xrange,yrange,xyFrac,Nx,Ny,FDorder,useSymmetry)
    % mesh = CreateMesh(xrange,yrange,Nx,Ny,FDorder,useSymmetry)
    % Creates uniform rectangular mesh from with ranges defined by Xrange
    % and Yrange, and Nx and Ny points in each direction. Finite difference
    % diference differentiation schemes with order FDorder are created.
    % If symmetry around the bottom y boundady, use usesymmetru=true. This
    % will move the bottom boundary points away from the axis.
    %
    % output : mesh is an object containing the mesh, differentiation
    % matrices and integration weights.
    
    tic;
    
    % Create mesh
    if (useSymmetry==true)
        % to use symmetry around the bottom y boundary, point at y = 0 is
        % removed
        dy =  yrange(2)/(Ny-.5);
        yrange(1) = yrange(1)+dy/2;       
    end
    
    x = linspace(xrange(1),xrange(2),Nx)';
    y = linspace(yrange(1),yrange(2),Ny)';
    dx = x(2)-x(1);
    dy = y(2)-y(1);

    [~,Nx1] = min(abs(x-xyFrac(1)));
    Nx2 = Nx;
    [~,Ny1] = min(abs(y-xyFrac(2)));
    Ny2 = Ny;
    
    x1  = x(    1:Nx1);
    x2  = x(Nx1+1:Nx2);
    xc  = mean([x(Nx1),x(Nx1+1)]);
    
    y1  = y(1:Ny1);
    y2  = y(1:Ny2);
    yc  = mean([y(Ny1),y(Ny1+1)]);
    
    %construct Mesh
    [Xl,Yl]  = meshgrid(x1,y1);
    [Xr,Yr]  = meshgrid(x2,y2);
        
    X = [Xl(:);Xr(:)];
    Y = [Yl(:);Yr(:)];
    
    Nxy = numel(X);
    
    %set up indexes
    indexes.leftDomain   = find(X(:)<xc);
    indexes.rightDomain  = find(X(:)>xc);
    indexes.bottomDomain = find(Y(:)<yc);
    indexes.topDomain    = find(Y(:)>yc);
    
    indexes.bottom    = find(Y(:)<min(Y(:))+dy/2);
    indexes.left      = find(X(:)<min(X(:))+dx/2) ;
    indexes.right     = find(X(:)>max(X(:))-dx/2) ;
    
    indexes.Lbottom   = find(abs(Y(:)-(yc-dy/2))<dy/2  & X(:)<xc+dx) ;
    indexes.Lright    = find(abs(X(:)-(xc+dx/2))<dx/2  & Y(:)>yc-dy) ;
    indexes.top       = find(Y(:)>max(Y(:))-dy/2& X(:)>xc  ) ;
    
    % 1D x derivative
    [Dxb_1D,D2xb_1D,~,~,~,~]                                         = Dmats_SBP(Nx2,dx,FDorder);
    [Dxt_1D,D2xt_1D,~,~,~,~]                                         = Dmats_SBP(Nx2-Nx1,dx,FDorder);
    % 1D y derivative
    [Dyl_1D,D2yl_1D,Dyl_1D_symm,Dyl_1D_asymm,D2yl_1D_symm,D2yl_1D_asymm] = Dmats_SBP(Ny1,dy,FDorder);
    [Dyr_1D,D2yr_1D,Dyr_1D_symm,Dyr_1D_asymm,D2yr_1D_symm,D2yr_1D_asymm] = Dmats_SBP(Ny2,dy,FDorder);
    
    %construct diff matrices for the domain
    Dxb          = sparse(kron(Dxb_1D ,speye(Ny1)));
    D2xb         = sparse(kron(D2xb_1D,speye(Ny1)));
    Dxt          = sparse(kron(Dxt_1D ,speye(Ny2-Ny1)));
    D2xt         = sparse(kron(D2xt_1D,speye(Ny2-Ny1)));
    
    Dx = sparse(Nxy,Nxy);
    D2x = sparse(Nxy,Nxy);
    p=indexes.bottomDomain;
        Dx(p,p) = Dxb;
        D2x(p,p) = D2xb;
    p=indexes.topDomain;
        Dx(p,p) = Dxt;
        D2x(p,p) = D2xt;
    
    %% Y Derivatives
    %construct diff matrices for the domain
    Dyl          = sparse(kron(speye(Nx1)    ,Dyl_1D ));
    Dyl_symm     = sparse(kron(speye(Nx1)    ,Dyl_1D_symm ));
    Dyl_asymm    = sparse(kron(speye(Nx1)    ,Dyl_1D_asymm ));

    D2yl         = sparse(kron(speye(Nx1)    ,D2yl_1D));
    D2yl_symm    = sparse(kron(speye(Nx1)    ,D2yl_1D_symm));
    D2yl_asymm   = sparse(kron(speye(Nx1)    ,D2yl_1D_asymm));
    
    Dyr          = sparse(kron(speye(Nx2-Nx1),Dyr_1D ));
    Dyr_symm      = sparse(kron(speye(Nx2-Nx1),Dyr_1D ));
    Dyr_asymm    = sparse(kron(speye(Nx2-Nx1),Dyr_1D_asymm ));
    
    D2yr         = sparse(kron(speye(Nx2-Nx1),D2yr_1D));
    D2yr_symm    = sparse(kron(speye(Nx2-Nx1),D2yr_1D_symm));
    D2yr_asymm   = sparse(kron(speye(Nx2-Nx1),D2yr_1D_asymm));
    
    Dy = sparse(Nxy,Nxy);
    D2y = sparse(Nxy,Nxy);
    p=indexes.leftDomain;
        Dy       (p,p)       = Dyl;
        Dy_symm  (p,p)       = Dyl_symm;
        Dy_asymm (p,p)       = Dyl_asymm;
        
        D2y      (p,p)       = D2yl;
        D2y_symm (p,p)       = D2yl_symm;
        D2y_asymm(p,p)       = D2yl_asymm;
        
    p=indexes.rightDomain;
    
        Dy       (p,p)       = Dyr;
        Dy_symm  (p,p)       = Dyr_symm;
        Dy_asymm (p,p)       = Dyr_asymm;
        
        D2y      (p,p)       = D2yr;
        D2y_symm (p,p)       = D2yr_symm;
        D2y_asymm(p,p)       = D2yr_asymm;
    
    
    % Pack outputs
    domains.n = 2;
    
    X1 = Xl;
    Y1 = Yl;
    X1(:,end+1) = X1(:,end)+dx;
    Y1(:,end+1) = Y1(:,end);
    indexes1  = find(X(:)<xc+dx & Y(:)<yc); 
    domains.X = {X1,Xr};
    domains.Y = {Y1,Yr};
    
    domains.indexes = {indexes1,indexes.rightDomain};
    mesh=struct('x',x,'y',y,'X',X,'Y',Y,'indexes',indexes,'domains',domains);
    for names = who('D*')' 
        eval([ 'mesh.' names{1} '=' names{1} ';' ]);
    end
    
    mesh.W = ones(size(X))*dx*dy;
    for p={indexes.top,indexes.left,indexes.right}
        mesh.W(p{1})=mesh.W(p{1})/2;
    end

    mesh.nGridPoints=numel(X);
    
    
disp(['    elapsed time - Creating uniform mesh: ' datestr(toc/24/3600, 'HH:MM:SS')]);



