function [x,u]=heattheta(xspan,tspan,nstep,theta,mu,...
               u0,g,f,varargin)
%HEATTHETA resuelve la ecuacion del calor con el
%  teta-metodo.
%  [X,U]=HEATTHETA(XSPAN,TSPAN,NSTEP,THETA,MU,U0,G,F)
%  resuelve la ecuacion del calor
%   D U/DT - MU D^2U/DX^2 = F en
%  (XSPAN(1),XSPAN(2)) X (TSPAN(1),TSPAN(2)) usando el
%  teta-metodo con condicion inicial U(X,0)=U0(X) y
%  condiciones de contorno de Dirichlet  U(X,T)=G(X,T)
%  para X=XSPAN(1) y X=XSPAN(2). MU es una constante
%  positiva, F, G y U0 funciones inline, NSTEP(1) es el
%  numero de intervalos de integracion espaciales,
%  NSTEP(2)+1 es el numero de intervalos de
%  integracion en tiempo.
h  = (xspan(2)-xspan(1))/nstep(1);
dt = (tspan(2)-tspan(1))/nstep(2);
N = nstep(1)+1;
 e =ones(N,1); D = spdiags([-e 2*e -e],[-1,0,1],N,N);
I = speye(N); A = I+mu*dt*theta*D/h^2;
An = I-mu*dt*(1-theta)*D/h^2;
 A(1,:) = 0;
 A(1,1) = 1; A(N,:) = 0; A(N,N) = 1;
 x =linspace(xspan(1),xspan(2),N);
 x = x'; fn = feval(f,x,tspan(1),varargin{:}); 
 un =feval(u0,x,varargin{:});
 [L,U]=lu(A); 
 for t = tspan(1)+dt:dt:tspan(2)
    fn1 = feval(f,x,t,varargin{:});
    rhs = An*un+dt*(theta*fn1+(1-theta)*fn);
    temp = feval(g,[xspan(1),xspan(2)],t,varargin{:});
    rhs([1,N]) = temp;
    u = L\rhs;
    u = U\u;
    fn = fn1;
    un = u;
	plot(x,u,'r')
%hold on
end
return
