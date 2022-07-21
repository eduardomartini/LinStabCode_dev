function [axs] = plotFlow(X,Y,vars,nr,nc,usedInd,varargin)
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
    %   usedInd (optional) : if the fields in vars do not contain all grid 
    %                       points (e.g., when masks are used), provide 
    %                       usedInd (typically mesh.usedInd to indicate
    %                       the grid points used.
    %   varargin (optional): extra arguments to contourf, used to create
    %                        plots. E.g. 'colorline','none'.
    % Outputs
    %   axs     : list with references to each of the subplots axis.
    

nVars   = size(vars,1);


if exist('nr');  if isempty(nr); nr = ceil(sqrt(nVars)); end; end
if exist('nc');  if isempty(nc); nc = nr               ; end; end

if ~( exist('usedInd') && ~isempty(usedInd)) 
    usedInd = 1:numel(X);
end

for i=1:nVars
    axs(i) = subplot(nr,nc,i);
    var = nan(size(X));
    
    var(usedInd) = vars{i,1};
    
    contour_varargin = {X,Y,var,varargin{:}};
    
    contourf(contour_varargin{:}); 
%     contourf(X,Y,var,'linecolor','none',varargin); 
    shading interp
%     pcolor(X,Y,var); shading interp
    colorbar;
%     axis equal tight
    title(vars{i,2});
end
