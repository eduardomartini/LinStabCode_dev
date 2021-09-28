function  DW = CreateDW(mesh,FDorder)
%CREATEDW Summary of this function goes here
%   Detailed explanation goes here
    if ~exist('FDorder')
        if isfield(mesh,'FDorder')
            FDorder=mesh.FDorder;
        else
            error('CreateDW : FD order not defined as an argument nor specified in the mesh object');
        end
    end
    
    N=mesh.ngp;
    [NX,NY] = size(mesh.X);
    
    Dx       = sparse(N,N);
    D2x      = sparse(N,N);
    Dy       = sparse(N,N);
    D2y      = sparse(N,N);
    Dy_symm  = sparse(N,N);
    D2y_symm = sparse(N,N);
    Dy_asymm = sparse(N,N);
    D2y_asymm= sparse(N,N);
    
    W       = zeros (NX,NY); 
    W_symm  = zeros (NX,NY); 
    W_asymm = zeros (NX,NY); 
    W      (mesh.usedInd)=1;
    W_symm (mesh.usedInd)=1;
    W_asymm(mesh.usedInd)=1;
    
    X = mesh.X(mesh.usedInd);
    Y = mesh.Y(mesh.usedInd);
   
    xList=sort(unique(X(:)));
    yList=sort(unique(Y(:)));
    
    dx = mesh.X(2,2)-mesh.X(1,1);
    dy = mesh.Y(2,2)-mesh.Y(1,1);
    
    ySym = dy/2;
    
    % Create X derivative Matrices
    
    for y = yList'
        % find all used points with a given x coordinate. 
        pos = find(Y==y);
        
        %guarantee that points are in increasing order
        [~,order] = sort(X(pos));
        pos=pos(order);
        
        %find holes 
        xl=X(pos);
        dx_line=(xl(2:end)-xl(1:end-1))/dx;
        domList = [0;find(dx_line>1.5);length(xl)];
        
        %create a FD scheme in each domain
        for idom = 1:length(domList)-1
           
            dom = domList(idom)+1:domList(idom+1);
            Nx  = length(dom);
            
            dom = domList(idom)+1:domList(idom+1);
            if Nx<2
                error('CreateDW : Cannot construct derivative with less than three points!')
            elseif Nx<FDorder+1
                FDorder_curr = 2;
            else
                FDorder_curr=FDorder;
            end
            
            [   Dx_1D     ,D2x_1D     , ...
                ~         ,     ~     , ...
                ~         ,     ~      ] = Dmats_SBP(Nx,dx,FDorder_curr);
        
            dom_pos = pos(dom);
            Dx      (dom_pos,dom_pos) = Dx_1D        ;
            D2x     (dom_pos,dom_pos) = D2x_1D       ;
          
            
            %trapeizodal rule
            
            meshpos=mesh.usedInd(dom_pos);
            W     (meshpos)          = W     (meshpos)*dx;
            W     (meshpos([1,end])) = W     (meshpos([1,end]))/2;
            
        end
    end
    
    % Create Y derivative Matrices

    for x = xList'
        % find all used points with a given x coordinate. 
        pos = find(X==x);
        
        %guarantee that points are in increasing order
        [~,order] = sort(Y(pos));
        pos=pos(order);
        
        %find holes 
        yl=Y(pos);
        dy_line=(yl(2:end)-yl(1:end-1))/dy;
        domList = [0;find(dy_line>1.5);length(yl)];
        
        %create a FD scheme in each domain
        for idom = 1:length(domList)-1
           
            dom = domList(idom)+1:domList(idom+1);
            Ny  = length(dom);

            if Ny<2
                error('CreateDW : Cannot construct derivative with less than three points!')
            elseif Ny<=FDorder+1
                FDorder_curr = 2;
            else
                FDorder_curr=FDorder;
            end
                
            [   Dy_1D      , D2y_1D     , ...
                Dy_1D_sym  , Dy_1D_asym , ...
                D2y_1D_sym , D2y_1D_asym ] = Dmats_SBP(Ny,dy,FDorder_curr);
        
            dom_pos = pos(dom);
            
            Dy      (dom_pos,dom_pos) = Dy_1D        ;
            D2y     (dom_pos,dom_pos) = D2y_1D       ;

            %trapeizodal rule
            meshpos=mesh.usedInd(dom_pos);
            W(meshpos)               = W(meshpos)*dy;
            W(meshpos([1,end]))      = W(meshpos([1,end]))/2;
            
            if (abs(yl(1)-ySym)<.1*dy) && mesh.symmetricBC
                Dy_symm  (dom_pos,dom_pos) = Dy_1D_sym    ;
                D2y_symm (dom_pos,dom_pos) = D2y_1D_sym   ;
                Dy_asymm (dom_pos,dom_pos) = Dy_1D_asym   ;
                D2y_asymm(dom_pos,dom_pos) = D2y_1D_asym  ;
                W_symm (meshpos)           = W_symm (meshpos)*dx;
                W_symm (meshpos(end))      = W_symm (meshpos(end))/2;
                W_asymm(meshpos)           = W_asymm(meshpos)*dx;
                W_asymm(meshpos( 1 ))      = W_asymm(meshpos( 1 )).*3/4;
                W_asymm(meshpos(end))      = W_asymm(meshpos(end))  ./2;
            else
                Dy_symm  (dom_pos,dom_pos) = Dy_1D        ;
                D2y_symm (dom_pos,dom_pos) = D2y_1D       ;
                Dy_asymm (dom_pos,dom_pos) = Dy_1D        ;
                D2y_asymm(dom_pos,dom_pos) = D2y_1D       ;
                W_symm (meshpos)           = W_symm (meshpos)*dx;
                W_symm (meshpos([1,end]))  = W_symm (meshpos([1,end]))/2;
                W_asymm(meshpos)           = W_asymm(meshpos)*dx;
                W_asymm(meshpos([1,end]))  = W_asymm(meshpos([1,end]))/2;
        end
    end
    
    
    
    DW.Dx  = Dx ; 
    DW.D2x = D2x; 

    DW.Dy  = Dy ;
    DW.D2y = D2y; 
    
    DW.Dxy = Dx*Dy ; 
    DW.Dyx = Dy*Dx ; 
    
    DW.W   = W ; 

    
    if mesh.symmetricBC
        DW.Dy_symm  = Dy_symm ;
        DW.D2y_symm = D2y_symm; 

        DW.Dxy_symm = Dx*Dy_symm ; 
        DW.Dyx_symm = Dy_symm*Dx ; 

        DW.Dy_asymm  = Dy_asymm ;
        DW.D2y_asymm = D2y_asymm; 

        DW.Dxy_asymm = Dx*Dy_asymm ; 
        DW.Dyx_asymm = Dy_asymm*Dx ; 

        DW.W_symm    = W_symm ; 
        DW.W_asymm   = W_asymm ; 

    else
        DW.Dy_symm  = Dy ;
        DW.D2y_symm = D2y; 

        DW.Dxy_symm = Dx*Dy; 
        DW.Dyx_symm = Dy*Dx ; 

        DW.Dy_asymm  = Dy;
        DW.D2y_asymm = D2y; 

        DW.Dxy_asymm = Dx*Dy; 
        DW.Dyx_asymm = Dy*Dx ; 
        
        DW.W_symm    = W; 
        DW.W_asymm   = W; 
        
    end
    

end

