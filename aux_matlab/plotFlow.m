function plotFlow(X,Y,vars,nr,nc)
    % plotFlow(X,Y,vars,nr,nc)
    %creates a subplot structure to make contour plots of flow quantities.
    % Inputs 
    %   X       : gridpoints x coordiantes
    %   Y       : gridpoints y coordiantes
    %   vars    : a N by 2 cell with each line containing  {F,label}, where
    %   F is a matrix with the same size of X and Y, to be ploted, and
    %   label is the title of the plot
    %   nr,nc   : number of rowns and columns of subplots. nr*nc >=
    %                   size(vars,1)
    

nVars   = size(vars,1);

if nargin==3
    nr      = ceil(sqrt(nVars));
    nc      = nr;
end

for i=1:nVars
    subplot(nr,nc,i)
    if size(vars{i,1})~=size(X)
        var     = reshape(vars{i,1},size(X));
    else
        var     = vars{i,1};
    end
    pcolor(X,Y,var); shading interp
    colorbar;
    axis equal tight
    title(vars{i,2});
end
