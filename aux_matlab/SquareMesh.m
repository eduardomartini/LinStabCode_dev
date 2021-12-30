function mesh = SquareMesh(xrange,yrange,Nx,Ny,FDorder,useSymmetry,periodic,alpha_filter)
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
    

    [X,Y]  = meshgrid(x,y);
    originalStructure=size(X);    
    
    % Prepare output
    mesh.X                 = X              ;
    mesh.Y                 = Y              ;
    mesh.ngp               = numel(X)       ;
    mesh.usedInd           = (1:numel(X))'  ;
    mesh.symmetricBC       = useSymmetry    ;
    mesh.FDorder           = FDorder        ;

    % Create dif matrices and int weights
    mesh.DW=CreateDW(mesh,periodic);


    %create low-pass filter
    if ~exist('alpha_filter')
        %no filtering
        mesh.alpga_filter='none'
        filter    = @(x) x;
        filter_ct = @(x) x;        
    else
        mesh.alpga_filter=alpha_filter;
        [filter,filter_ct]  = GetFilter(mesh,alpha_filter);
    end
    mesh.Filters.filter    =  filter    ;
    mesh.Filters.filter_ct =  filter_ct ;
    
    % get boundary indexes
    indexes.bottom    = find(Y(mesh.usedInd)==min(Y(:)));
    indexes.top       = find(Y(mesh.usedInd)==max(Y(:)));
    indexes.left      = find(X(mesh.usedInd)==min(X(:)));
    indexes.right     = find(X(mesh.usedInd)==max(X(:)));
    mesh.indexes=indexes;