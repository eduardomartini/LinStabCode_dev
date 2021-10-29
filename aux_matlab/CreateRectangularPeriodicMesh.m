function mesh = CreateRectangularPeriodicMesh(xrange,yrange,Nx,Ny,FDorder,useSymmetry)
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
    
    x = linspace(xrange(1),xrange(2),Nx+1)';
    x = x(1:end-1);
    y = linspace(yrange(1),yrange(2),Ny)';
    
    if Nx == 1
        dx=1
    else
        dx = x(2)-x(1);
    end
    dy = y(2)-y(1);

    [X,Y]  = meshgrid(x,y);
    
    

    % 1D x derivative
    
    %Create matri to impose periodicity
    if Nx == 1 
        Dx_1D  = 0;
        D2x_1D = 0;
    else
        B = sparse(Nx,Nx+FDorder*2);
        for i=1:Nx
            B(i,i+FDorder)=1;
        end
        for i=1:FDorder
            B(end+(i-1)-(FDorder-1),i)=1;
            B(i,end+(i-1)-(FDorder-1))=1;

        end   
        %%
        [Dx_1D,D2x_1D,~,~,~,~]                                         = Dmats_SBP(Nx+FDorder*6,dx,FDorder);
        Dx_1D_periodic = Dx_1D ( (1:Nx)+FDorder*2,(1:Nx+FDorder*2)+FDorder);
        Dx_2D_periodic = D2x_1D( (1:Nx)+FDorder*2,(1:Nx+FDorder*2)+FDorder);


        Dx_1D = Dx_1D_periodic *B';
        D2x_1D= Dx_2D_periodic *B';
    end
    
    % 1D y derivative
    [Dy_1D,D2y_1D,Dy_1D_symm,Dy_1D_asymm,D2y_1D_symm,D2y_1D_asymm] = Dmats_SBP(Ny,dy,FDorder);
    
    %construct diff matrices for the domain
    Dx          = sparse(kron(Dx_1D,speye((Ny))));
    D2x         = sparse(kron(D2x_1D,speye(Ny)));

    Dy          = sparse(kron(speye((Nx)),Dy_1D));
    Dy_symm     = sparse(kron(speye((Nx)),Dy_1D_symm));
    Dy_asymm    = sparse(kron(speye((Nx)),Dy_1D_asymm));
    D2y         = sparse(kron(speye(Nx),D2y_1D));
    D2y_symm    = sparse(kron(speye((Nx)),D2y_1D_symm));
    D2y_asymm   = sparse(kron(speye((Nx)),D2y_1D_asymm));

    Dxy         = Dy        *Dx      ;
    Dyx         = Dx        *Dy      ;
    Dxy_symm    = Dy_symm   *Dx      ;
    Dyx_symm    = Dx        *Dy_symm ;
    Dxy_asymm   = Dy_asymm  *Dx      ;
    Dyx_asymm   = Dx        *Dy_asymm;
    
    indexes.bottom    = find(Y(:)<min(Y(:))+dy/2);
    indexes.top       = find(Y(:)>max(Y(:))-dy/2);
    indexes.left      = find(X(:)<min(X(:))+dx/2);
    indexes.right     = find(X(:)>max(X(:))-dx/2);
   
    
    
    % Pack outputs
    mesh=struct('X',X,'Y',Y);
    for names = who('D*')' 
        eval([ 'mesh.DW.' names{1} '=' names{1} ';' ]);
    end
    mesh.usedInd = (1:numel(X))';
    mesh.ngp = numel(X);

    mesh.DW.W = ones(size(X))*dx*dy;
    for p={indexes.top,indexes.bottom}
        mesh.DW.W(p{1})=mesh.DW.W(p{1})/2;
    end
    mesh.nGridPoints=numel(X);
    
disp(['    elapsed time - Creating uniform mesh: ' datestr(toc/24/3600, 'HH:MM:SS')]);



