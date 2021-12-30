function [maskedMesh,edges] = MeshMask(mesh,mask)
%MESHMASK Summary of this function goes here
%   Detailed explanation goes here
    maskedMesh                   = mesh;
    maskedMesh.usedInd           = find(~mask);
    maskedMesh.ngp               = numel(maskedMesh.usedInd);
    maskedMesh.DW                = CreateDW(maskedMesh);
    
    %find edges
    maskedge=mask;
    maskedge(2:end  ,:) =  maskedge(2:end  ,:)+mask(1:end-1,:);
    maskedge(1:end-1,:) =  maskedge(1:end-1,:)+mask(2:end  ,:);
    maskedge(:,2:end  ) =  maskedge(:,2:end  )+mask(:,1:end-1);
    maskedge(:,1:end-1) =  maskedge(:,1:end-1)+mask(:,2:end  );

    edges = find(maskedge(maskedMesh.usedInd)>0);
    
    %% Update filter
    if mesh.alpga_filter=='none'
        %no filtering
        filter    = @(x) x;
        filter_ct = @(x) x;        
    else
        [filter,filter_ct]  = GetFilter(mesh,alpha_filter);
    end
    mesh.Filters.filter    =  filter    ;
    mesh.Filters.filter_ct =  filter_ct ;

end

