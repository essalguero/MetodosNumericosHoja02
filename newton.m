function [zero,res,niter]=newton(fun,dfun,x0,tol,...
                                 nmax,varargin)
%NEWTON Trova uno zero di una funzione.
%   ZERO=NEWTON(FUN,DFUN,X0,TOL,NMAX) approssima lo
%   zero ZERO della funzione definita nella function
%   FUN, continua e derivabile, usando il metodo di
%   Newton e partendo da X0. Se la ricerca
%   dello zero fallisce, il programma restituisce un
%   messaggio d'errore. FUN e DFUN possono essere
%   inline function, anonymous function o function
%   definite in M-file.
%   ZERO=NEWTON(FUN,DFUN,X0,TOL,NMAX,P1,P2,...) passa
%   i parametri P1,P2,... alle funzioni
%   FUN(X,P1,P2,...) e DFUN(X,P1,P2,...).
%   [ZERO,RES,NITER]= NEWTON(FUN,...) restituisce il
%   valore del residuo RES in ZERO ed il numero di
%   iterazioni NITER necessario per calcolare ZERO.

x = x0;
fx = fun(x,varargin{:});
dfx = dfun(x,varargin{:});
niter = 0; diff = tol+1;
while diff >= tol && niter < nmax
   niter = niter + 1;      diff = - fx/dfx;
   x = x + diff;           diff = abs(diff);
   fx = fun(x,varargin{:});
   dfx = dfun(x,varargin{:});
end
if (niter==nmax && diff > tol)
  fprintf(['Newton si e'' arrestato senza aver ',...
   'soddisfatto l''accuratezza richiesta, avendo\n',...
   'raggiunto il massimo numero di iterazioni\n']);
end
zero = x; res = fx;
return
