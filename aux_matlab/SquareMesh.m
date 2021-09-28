function mesh = SquareMesh(xrange,yrange,Nx,Ny,FDorder,useSymmetry)
    % mesh = CreateMesh(xrange,yrange,Nx,Ny,FDorder,useSymmetry)
    % Creates uniform rectangular mesh from with ranges defined by Xrange
    % and Yrange, and Nx and Ny points in each direction. Finite difference
    % diference differentiation schemes with order FDorder are created.
    % If symmetry around the bottom y boundady, use usesymmetru=true. This
    % will move the bottom boundary points away from the axis.
    %
    % output : mesh is an object containing the mesh, differentiation
    % matrices and integration weights.
        
    % Create mesh
    if (useSymmetry==true && yrange(1)==0)
        % to use symmetry around the bottom y boundary, point at y = 0 is
        % removed
        dy =  yrange(2)/(Ny-.5);
        yrange(1) = yrange(1)+dy/2;       
    end
    
    x = linspace(xrange(1),xrange(2),Nx)';
    y = linspace(yrange(1),yrange(2),Ny)';
    
    dx = x(2)-x(1);
    dy = y(2)-y(1);

    [X,Y]  = meshgrid(x,y);
    originalStructure=size(X);    
    
    % Prepare output
    mesh.X                 = X          ;
    mesh.Y                 = Y          ;
    mesh.ngp               = numel(X)   ;
    mesh.usedInd           = 1:numel(X) ;
    mesh.symmetricBC       = useSymmetry;
    mesh.FDorder           = FDorder;

    mesh.DW=CreateDW(mesh);


