function [u,x,y,error]=poissonfd(a,c,b,d,nx,ny,fun,...
                       bound,uex,varargin)
%POISSONFD resolve o problema de Poisson bidimensional
%  [U,X,Y]=POISSONFD(A,C,B,D,NX,NY,FUN,BOUND) resolve
%  pelo esquema de diferenas finitas de cinco pontos,
%  o problema -LAPL(U) = FUN no rectngulo (A,C)X(B,D)
%  com condies na fronteira de Dirichlet U(X,Y)=
%  BOUND(X,Y) para qualquer (X,Y) na fronteira do
%  rectngulo.
%
%  [U,X,Y,ERROR]=POISSONFD(A,C,B,D,NX,NY,FUN,BOUND,UEX)
%  calcula tambm o mximo do erro nodal, ERROR, em
%  relao  soluo exacta UEX. FUN, BOUND e UEX
%  podem ser funes inline.
if nargin == 8
    uex = inline('0','x','y');
end
nx=nx+1;   ny=ny+1;  hx=(b-a)/nx; hy=(d-c)/ny;
nx1=nx+1;  hx2=hx^2; hy2=hy^2;
kii=2/hx2+2/hy2;     kix=-1/hx2;  kiy=-1/hy2;
dim=(nx+1)*(ny+1);   K=speye(dim,dim);
rhs=zeros(dim,1);
y = c;
for m = 2:ny
    x = a; y = y + hy;
    for n = 2:nx
        i = n+(m-1)*(nx+1);
        x = x + hx;
        rhs(i) = feval(fun,x,y,varargin{:});
        K(i,i) = kii;         K(i,i-1) = kix;
        K(i,i+1) = kix;       K(i,i+nx1) = kiy;
        K(i,i-nx1) = kiy;
    end
end
rhs1 = zeros(dim,1);
x = [a:hx:b];
rhs1(1:nx1) = feval(bound,x,c,varargin{:});
rhs1(dim-nx:dim) = feval(bound,x,d,varargin{:});
y = [c:hy:d];
rhs1(1:nx1:dim-nx) = feval(bound,a,y,varargin{:});
rhs1(nx1:nx1:dim) = feval(bound,b,y,varargin{:});
rhs = rhs - K*rhs1;
nbound = [[1:nx1],[dim-nx:dim],...
          [1:nx1:dim-nx],[nx1:nx1:dim]];
ninternal = setdiff([1:dim],nbound);
K = K(ninternal,ninternal);
rhs = rhs(ninternal);
utemp = K\rhs;
uh = rhs1;
uh (ninternal) = utemp;
k = 1; y = c;
for j = 1:ny+1
    x = a;
    for i = 1:nx1
        u(i,j) = uh(k);
        k = k + 1;
        ue(i,j) = feval(uex,x,y,varargin{:});
        x = x + hx;
    end
    y = y + hy;
end
x = [a:hx:b];
y = [c:hy:d];
if nargout == 4
  if nargin == 8
     warning('Soluo exacta no disponvel');
     error = [ ];
  else
     error = max(max(abs(u-ue)))/max(max(abs(ue)));
  end
end
return
