function  [B,F,laKMM,errorr]=slnmf(k,lambda,iter,M)
  %%%%The model----------------
% min{B,F}w1*||M-BF||^2+lambda*Tr(Z'LYZ), 
% s.t. B>=0,F>=0, Z'Z=I.
% You only need to provide the above five inputs.
% Notation:
% M ... (n x n) cell-cell network 
%       n  ... number of cells            
% d ... number of features for clustering in NMF(the number of cluster)
% gamma ... Regularization parameter
%
% iter ... The maximum number of iterations
%%%
% Coder Wenming Wu Email: wenmingwu55 at 163.com
%
%   version  --December/2020

    errorr=zeros(iter,1);
    [~,n]=size(M);
    l=k;
    %initialization for variables B,F 
	[U,V,D]=svds(M,k);
	%normalization on it 
    B=abs(U*sqrt(V));
    F=abs(sqrt(V)*D');
    P=zeros(n,n);
    Pk=zeros(k,k);
    L=[P,B;B',Pk]; 
	
    Ds=diag(sum(L,2));
	Ls=Ds-L;
	
    [a,aa]=eig(Ls);
	a=[a;diag(aa)';abs(diag(aa))']';
    a=sortrows(a,n+k+2);
	Ft=a(2:l+1,1:n+k)';
	V=zeros(n,l);
	
    for kkk=1:n
    for oo=1:k
	V(kkk,oo)=norm(Ft(kkk,:)-Ft(n+oo,:));
    end
    end
    D1=1;
    D2=1;
    %%%%%%%===========Update variables B,F by iteration================
for o=1:iter
    U1=D1*U;
	V1=D2*V;
	dist=sqdist(U1',V1');
	tmp=zeros(n,k);
    for i=1:n
        B(i,:)=B(i,:).*((2*M(i,:)*F')./(2*B(i,:)*(F*F')+lambda*V(i,:))); 
    end
    F=F.*((B'*M)./(B'*B*F+eps));

    L=[P,B;B',Pk];
	Ds=diag(sum(L,2));
	Ls=Ds-L;
	[a,aa]=eig(Ls);
	a=[a;diag(aa)';abs(diag(aa))']';
    a=sortrows(a,n+k+2);
	Ft=a(1:l,1:n+k)';
	V=zeros(n,l);
    for kkk=1:n
        for oo=1:k
            V(kkk,oo)=norm(Ft(kkk,:)-Ft(n+oo,:));
        end
    end
    
    errorr(o,1) = norm(M-B*F,'fro');
    if norm(M-B*F,'fro')<0.01
        break;
    else 
        B=B;
        F=F;
    end
        [clusternum, laKMM] = struG2la(B);
        [Ct,~]=kmeans(B,l,'Replicates',100);

end
end
