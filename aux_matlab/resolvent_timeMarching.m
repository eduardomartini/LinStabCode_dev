function [S,U,V] = resolvent(L0,W,invW,omega,nEig,B,C,filters,dt,timeMarchMethod)
% [gain,U_out,V_in] = ResAnalysis(LHS,RHS,F,omega,noEigs,B,C)
% Performs classical resolvent analysis using "eigs" function of the
% linear operator H=(L0-i*omega)^-1 using input and output matrices
% given by B and C. The adjoint is obtained with respect to the norm W.
% Inputs :
%     noEigs : number of eigen values desired.
% Outputs :
%     S : singular values (sqrt. of gains)
%     U : optimal outputs/responses
%     V : optimal inputs/forcing

% reduce to generalized to standard EVP
nDOF    = size(L0,1);

disp(['    computing resolvent analysis for omega=' num2str(omega)]);

if size(W,1)==1 || size(W,2)==1
    invW    = 1./W;
    W       = diag(sparse(  W ));
    invW    = diag(sparse(invW));
else
    if ~exist('invW') || prod(size(W) == size(invW))==0
        error('W is a matrix and its inverse was not provided or does not has the same size as W.');
    end
end

% LU-decomposition of L0-sigma*I
tic
if strcmp(timeMarchMethod,'backwardsEuler')
    Al = speye(size(L0))-dt*L0;
    Ar = speye(size(L0));
else
    error('Time marching method invalid')
end
[iLL,iUU]    = ilu(Al,struct('type','ilutp','droptol',1e-5));
time = toc;
disp(['    elapsed time - iLU-decomposition of L-sigma*I: ' datestr(time/24/3600, 'HH:MM:SS')]);

% Create spatial filter function from filters in the mesh object
if exist('filters')
    filter    = filters.filter;
    filter_ct = filters.filter_ct;
    FILTER    = @(x) reshape( filter   ( reshape(x,[],5)),[],1) ; 
    FILTER_ct = @(x) reshape( filter_ct( reshape(x,[],5)),[],1) ; 
else
    FILTER       = @(x)x;
    FILTER_ct    = @(x)x;
end

% Function handle for resolvent operator H and Hermitian H*H
invAl    = @(x) gmresnomsg(Al ,x,[],1e-6,100,iLL,iUU);
invAl_ct = @(x) gmresnomsg(Al',x,[],1e-6,100,iUU',iLL');


q0 = zeros(nDOF,1);
Dir = @(x) TimeMarch(q0,dt*x,invAl,Ar,1e-6/dt);
Adj = @(x) TimeMarch(q0,dt*x,invAl_ct,Ar',1e-6/dt);
HtrH        = @(v) invW*Adj(W*Dir(v));

% 'eigs' parameters
opts.tol    = eps;
opts.disp   = 2;
opts.issym  = false;
opts.isreal = false;

% Input via eigendecomposition of H*H
tic
[V, OMEGA] 	= eigs(HtrH, nDOF, nEig, 'lm', opts);
S           = sqrt(diag(real(OMEGA))); % singular values of H given by sqrt. of eigenvalues of H*H, should be real (imag. part is O(eps))

time = toc;
disp(['    elapsed time - SVD via ''eigs'' using matrix-free Arnoldi: ' datestr(time/24/3600, 'HH:MM:SS')]);

% Normalization and outputs
U = zeros(size(V));
for i = 1:nEig
    % normalization in inner product norm
    V(:,i)  = V(:,i)/sqrt(V(:,i)'*W*V(:,i));
%     U(:,i) =  H(V(:,i))/S(i); %C*(qq*(UU\(LL\(pp*(rr\(B*(V(:,i))))))))/S(i);
end
