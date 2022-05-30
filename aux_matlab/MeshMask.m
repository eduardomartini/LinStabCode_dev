function maskedMesh = MeshMask(mesh,mask)
%MESHMASK Summary of this function goes here
%   Detailed explanation goes here
    maskedMesh                   = mesh;
    maskedMesh.usedInd           = find(~mask);
    maskedMesh.ngp               = numel(maskedMesh.usedInd);
    maskedMesh                   = AddDW2mesh(maskedMesh);
    
    %find edges
    neighbours = ones(3,3);
    maskedge = conv2(mask,neighbours,'same');
    
  
    
    % update boundary indexes 
    for bondaries = fields(mesh.idx)'
        domain = zeros(size(mesh.X));
        ids = mesh.idx.(bondaries{1});
        domain(ids)=1;
        maskedMesh.idx.(bondaries{1}) = find(domain(maskedMesh.usedInd)==1);
    end
    
    %add mask edge
    edges = find(maskedge(maskedMesh.usedInd)>0);
    maskedMesh.idx.mi = edges;
    
    %% Update filter
    if mesh.alpga_filter=='none'
        %no filtering
        filter    = @(x) x;
        filter_ct = @(x) x;        
    else
        [filter,filter_ct]  = GetFilter(maskedMesh,maskedMesh.alpha_filter);
    end
    maskedMesh.filters.filter    =  filter    ;
    maskedMesh.filters.filter_ct =  filter_ct ;

end

