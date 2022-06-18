n=100;
A = diag((1:n).^1.5);
for i=1:n-1
    A(i,i+1) = 1;
end
fun = @(x) A*x;

neigs   = 3;
nkirlov = 21;
mkirlov = 1 ; 
pkirlov = 2 ;
tol     = 1e-8;

v = randn(n,1);
for i=1:n
    [v,~] = qr(A*v,0);
end
%%
w=A*v
[norm(w-v*v'*w),v'*w,((n-1)/n)^n]
%%
% tol     = 1e-3;
% [v,l] = reigs(fun,n,neigs,nkirlov,mkirlov,pkirlov,tol);
[v,l,r]  = eigs_blockshur(fun,n,neigs,nkirlov,mkirlov,pkirlov,tol);
[~,~,r2] = eigs_blockshur(fun,n,neigs,nkirlov,10,pkirlov,tol);
% % [l,diag(eig(g(v'*A*v)]
%%
[V,L] = eig(A);
[L,order]=sort(diag(L),'descend');
V = V(:,order);

e = v(:,1:neigs) - V(:,1:neigs);
err = diag(e'*e);
el = L(1:neigs)-l(1:neigs);
disp([el,err])

ncalls = size(r,2);
x=1:ncalls;

ncalls = size(r2,2);
x2=1:ncalls;
semilogy(x,r,'k',x2,r2,'bo',x,x*0+tol,'r')
% 
% vv=v(:,1:neigs)
% w = A*vv;
% H = vv'*w;
% fk = w-vv*H;
% % fk'*fk
