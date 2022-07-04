function [L,Lw,R,idx,L0,R0,R1,R2] = GetSpatialLinProblem(mesh,BF,model,floquetExp)
    % Constructs spatial linear operator
    % (L + omega Lw] q' = alpha R q',
    % where q' = [q; alpha q]
    if ~exist('floquetExp','var'); floquetExp=0 ;end

    tic; 
    disp('Computing spatial stability operator')
    if ~strcmp(model,'2D')
        error('Currently spatial stability is only valid for 2d analysis')
    end
    alpha_0  =  0;
    alpha_p1 =  1;
    alpha_m1 = -1;
    
    
    % RHS : R0  + alpha R1 + alpha^2 R2
    % Using solutions for RHS with alpha 0,-1 and +1 to obtain Ri
    [A_0 ,idx] = GetLinProblem(mesh,BF,model,alpha_0 ,floquetExp);
    [A_p1,~]   = GetLinProblem(mesh,BF,model,alpha_p1,floquetExp);
    [A_m1,~]   = GetLinProblem(mesh,BF,model,alpha_m1,floquetExp);
    n = size(A_0,1);
    Z = sparse(n,n);
    I = speye(n);
    
    L0 = speye(n)*-1i;
    R0 =  A_0;
    R1 = (A_p1 - A_m1)/2;
    R2 = (A_p1 - R0 - R1);
    
    % From -iw q = (R0 + alpha R1 + alpha^2 R2), construct
    % [ 0              I   ]    [         q ]          [  I      0  ]  [         q ] 
    % [ \omega I-R0        ]    [ alpha   q ]  = alpha [  R1     R2 ]  [ alpha   q ]

    L  = [ Z , I ;-R0 , Z ] ;
    Lw = [ Z , Z ; I  , Z ]*-1i ;
    
    R = [I,Z;R1,R2];
    
    %update idx to include extended space
    for f = fields(idx)
        ids = idx.(f{1});
        idx.(f{1}) = [ids,ids+n];
    end
    
    disp(['    elapsed time - Spatial stability matrices:',datestr(toc/24/3600, 'HH:MM:SS')]);

    