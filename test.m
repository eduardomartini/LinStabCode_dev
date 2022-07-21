addpath('aux_matlab');
opts.tol    = 1e-8;%sqrt(eps);
opts.disp   = 2;
opts.issym  = false;
opts.isreal = false;
opts.p      = 150;


n=5e3;
L1 = spdiags((1:n)'+1i,0,n,n);
L1 = sprand(n,n,0.1);

lambdatar=-1;
nEig=50;

% 

tic
[V,lambda] = eigs(L1,nEig,lambdatar,opts);
lambda=diag(lambda);

% [invA_fun,invA_T_fun] = GetInverseFunction(L1-speye(size(L1))*lambdatar);
% [V2,lambda2] = eigs(invA_fun,size(L1,1),nEig,lambdatar,opts); 
% lambda2=diag(lambda2);

time_eigsMatlab = toc;
tic
n = 500;
m = 104;        
p = 4;

% opts.type='builtin';
clear opts_invf
opts_invf.type='ilu';

%%
tic
[op,~] = GetInverseFunction(L1-speye(size(L1))*lambdatar,opts_invf);
[V4,ritz4] = block_eigs2(op, size(L1,1), nEig,n,m,p,opts.tol);
time_eigsblock3 = toc;
lambda4 = 1./(ritz4)+lambdatar;


%%
figure
plot(real(lambda),imag(lambda),'g*',real(lambda4),imag(lambda4),'sk')
disp([time_eigsMatlab,time_eigsblock3])


%%
% 
% opts_inf_f.type='ilu';
% opts_inf_f.verbose=true;
% opts_inf_f.tol=1e-8;
% [~,op] = GetInverseFunction(L1-speye(size(L1))*lambdatar,opts_inf_f);
% 
% opts_inf_f.type='lu';
% [~,op2] = GetInverseFunction(L1-speye(size(L1))*lambdatar,opts_inf_f);
% opts_inf_f.type='builtin';
% [~,op3] = GetInverseFunction(L1-speye(size(L1))*lambdatar,opts_inf_f);
% 
% V=randn(size(L1,1),1);
% norm(op(V)-op2(V))
% norm(op3(V)-op2(V))