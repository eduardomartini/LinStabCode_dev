function [Wen,invWen] = GetXuEnergyNorm(mesh,BF,model)

if model == 'axy'
    intWeight = mesh.DW.W .*mesh.Y; %int () r dr dz
elseif  model == '2D'
    intWeight = mesh.DW.W ; %int () dy dz
end    

p=mesh.usedInd;
Wen  = ...
    0.5*[
    (BF.T(p) ./(BF.RHO(p)*BF.kappa*BF.Ma^2)).* intWeight(p)
     BF.RHO(p).*intWeight(p)
     BF.RHO(p).*intWeight(p)
     BF.RHO(p).*intWeight(p)
    (BF.RHO(p)./(BF.kappa*(BF.kappa-1)*BF.T(p)*BF.Ma.^2)).*intWeight(p)
    ];

nDOFs = mesh.ngp*5;
invWen  = spdiags(1./Wen,0,nDOFs,nDOFs);
Wen     = spdiags(Wen   ,0,nDOFs,nDOFs);
% invF = 1./F;