function BF = ReadStepBaseflow(mesh,BF)
    
    X = mesh.X;
    Y = mesh.Y;
    
    BFdata = load('BaseFlows/Step_BF')  ; %reads 'data' from mat file
    
    BF.RHO      = X*0+1;    % Baseflow Density 
    BF.U        = X*0;      % Baseflow U velocity
    BF.V        = X*0;      % Baseflow V velocity
    BF.W        = X*0;      % Baseflow W velocity
    BF.T        = X*0+1;    % Baseflow Temperature
    BF.MU       = X*0+1;    % Baseflow viscosity 
    BF.cv       = 1/(BF.kappa*(BF.kappa-1)*BF.Ma^2); % ideal gas ???
    BF.c1       = (BF.kappa-1)*BF.Re*BF.Pr*BF.Ma^2;  % ideal gas ???
    BF.c2       = BF.kappa*BF.Ma^2;
    
    % derivatives of the viscosity can be provided here, or overwriten
    % latter using, e.g., the sutherland_air function.
    BF.dMUdT    = X*0;      % Baseflow viscosity variation with temperature 
    BF.d2MUdT2  = X*0;      % Baseflow viscosity 2nd derivative wrt temperature
    BF.kappa    = 1.4;      % ideal gas adiabatic coefficient

    % Interpolate W and V velocities from database.
    BF.W(mesh.usedInd) = griddata(BFdata.X,BFdata.Y,BFdata.U,X(mesh.usedInd),Y(mesh.usedInd)); 
    BF.U(mesh.usedInd) = griddata(BFdata.X,BFdata.Y,BFdata.V,X(mesh.usedInd),Y(mesh.usedInd)); 

    nanMask = nan(size(X));
    nanMask(mesh.usedInd)=1;
    
    BF.RHO      = BF.RHO.*nanMask;    % Baseflow Density 
    BF.U        = BF.U  .*nanMask;      % Baseflow U velocity
    BF.V        = BF.V  .*nanMask;      % Baseflow V velocity
    BF.W        = BF.W  .*nanMask;      % Baseflow W velocity
    BF.T        = BF.T  .*nanMask;    % Baseflow Temperature
    BF.MU       = BF.MU .*nanMask;    % Baseflow viscosity 
    
