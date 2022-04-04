function [LHS,RHS] = SetBoundaryConditions(LHS,RHS,indexes,borders,variables)   
tic
% inlet: zero value
index_set = [] ;
for b=borders
    for v=variables 
        index_str = [b,'i_' v];
        index_set = [index_set,indexes.(index_str)];
    end
end

index_set= unique(index_set(:));

LHS(index_set, :)   = 0;
RHS(index_set, :)   = 0;

LHS(index_set, index_set)   = eye(length(index_set));
RHS(index_set, index_set)   = eye(length(index_set));
    
time    = toc;
disp(['    elapsed time - Boundary conditions: ' datestr(time/24/3600, 'HH:MM:SS')]);

