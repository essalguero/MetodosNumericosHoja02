function [t,u]=feuler(odefun,tspan,y0,Nh,varargin)
%FEULER  Risolve equazioni differenziali
%   usando il metodo di Eulero in avanti.
%   [T,Y] = FEULER(ODEFUN,TSPAN,Y0,NH) con
%   TSPAN = [T0,TF] integra il sistema di equazioni
%   differenziali y' = f(t,y) dal tempo T0 a TF con
%   condizione iniziale Y0 usando il metodo di Eulero
%   in avanti su una griglia equispaziata di NH
%   intervalli.
%   La funzione ODEFUN(T,Y) deve ritornare un vettore
%   contenente f(t,y), della stessa dimensione di y.
%   Ogni riga del vettore soluzione Y corrisponde ad
%   un istante temporale del vettore colonna T.
%   [T,Y] = FEULER(ODEFUN,TSPAN,Y0,NH,P1,P2,...) passa
%   i parametri addizionali P1,P2,... alla funzione
%   ODEFUN come ODEFUN(T,Y,P1,P2...).
h=(tspan(2)-tspan(1))/Nh;
y=y0(:); % genera sempre un vettore colonna
w=y; u=y.';
tt=linspace(tspan(1),tspan(2),Nh+1);
for t = tt(1:end-1)
 w=w+h*odefun(t,w,varargin{:});
 u = [u; w.'];
end
t=tt';
return
