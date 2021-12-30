function BF = ReadBaseflow(mesh,Re,m,Ma,Pr)
    
    X = mesh.X;
    Y = mesh.Y;
    
    load Paraboloid_BF; %reads 'data' from mat file
    
    Xnek = data(:,:,1);
    Ynek = data(:,:,2);
    Unek = data(:,:,3);
    Vnek = data(:,:,4);
    clear data;

    BF.Ma = Ma;
    BF.Pr = Pr;
    
    BF.RHO      = X*0+1;
    BF.U        = X*0;
    BF.V        = X*0;
    BF.W        = X*0;

    BF.W(:) = griddata(Xnek(:),Ynek(:),Unek(:),X(:),Y(:));
    BF.U(:) = griddata(Xnek(:),Ynek(:),Vnek(:),X(:),Y(:));

    BF.T        = X*0+1;
    BF.MU       = X*0+1;
    BF.kappa    = 1.4;
    BF.cv       = 1/(BF.kappa*(BF.kappa-1)*Ma^2);
    BF.c1       = (BF.kappa-1)*Re*Pr*Ma^2;
    BF.c2       = BF.kappa*Ma^2;


    BF.dMUdT = X*0;
    BF.dmudT = X*0;

    BF.d2MUdT2 = X*0;
    BF.d2mudT2 = X*0;

    BF.m = m;
    BF.Re = Re;