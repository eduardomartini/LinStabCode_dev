function x = quiet_gmres(varargin)
    [x,~] = gmres(varargin{:});
end
