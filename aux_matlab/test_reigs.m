n=100;
A = diag((1:n).^1.5);
for i=1:n-1
    A(i,i+1) = 1;
end
fun = @(x) A*x;

neigs   = 10;
nkirlov = 11;
mkirlov = 1 ; 
pkirlov = 2 ;
tol     = 1e-8;

%%
[v,l,r]  = eigs_blockshur(fun,n,neigs,nkirlov*pkirlov,mkirlov,1,tol);
[~,~,r2] = eigs_blockshur(fun,n,neigs,nkirlov,mkirlov,pkirlov,tol);

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
subplot(211)
    semilogy(x,r,'k',x2,r2,':b',x,x*0+tol,'r')
    xlabel('Iterations')
subplot(212)
    semilogy(x,r,'k',x2*pkirlov,r2,':b',x,x*0+tol,'r')
    xlabel('Function calls')
% 
% vv=v(:,1:neigs)
% w = A*vv;
% H = vv'*w;
% fk = w-vv*H;
% % fk'*fk
