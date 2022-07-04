function [L0,index_set] = SetBoundaryConditions(L0,idx,borders,variables,spatialprob)   
    % L0,index_set] = SetBoundaryConditions(L0,idx,borders,variables)   
    % Imposes Dirichllet boundary conditions on L0 at the selected borders
    % for the choosen variables
    % Inputs :
    %   L0          : the linear operator to be changed
    %   idx         : a structure containing the indexes of each border for
    %                 each variable
    %   borders     : string array listing the borders to which the b.c. is 
    %                 to be applied. For standanrt mesh, 'lrbt' indicate the 
    %                 left, right, bottom and top boundaries
    %   variables   : string array listing the variables to which the b.c. 
    %                 are to te applied. Options,  ruvwT', for density (rho),
    %                 the u,v,w velocities and temperature.
    %  spatialprob  : if true, matrix is assumed to correspond to a spatial 
    %                 stability problem.

if ~exist('spatialprob','var') ; spatialprob=false; end

tic
% inlet: zero value
index_set = [] ;
for b=borders
    for v=variables 
        if v=='r';  v='rho'; end  % if density, expand the string
        index_str = [b,'i_' v];
        index_set = [index_set;idx.(index_str)(:)];
    end
end

index_set= unique(index_set(:));

if spatialprob % if spatial problem, add BC to q and to alpha q
    n                   = size(L0);
    index_set           = [index_set(:),index_set(:)+n];
    L0(index_set, :)    = 0;
end

L0(index_set, :) = 0;
L0(:, index_set) = 0;
L0(index_set, index_set)   = eye(length(index_set));
time    = toc;

disp(['    elapsed time - Boundary conditions: ' datestr(time/24/3600, 'HH:MM:SS')]);

