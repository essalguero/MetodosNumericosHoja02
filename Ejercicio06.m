clear all
clc

format long

T = 30; L = 2; alpha = 0.004

theta = 0; %Euler Explicito

f = @(x, t) 2 - x + exp(-2 * t)
g = @(x, t) (2 / (1 + (t - (1 / 2))^2)) - x.*(2 - x)

u0 = inline('2.*(1 - x).^2.*(heaviside(x-1/2)-heaviside(x-3/2))', 'x')

xspan = [0, L]
tspan = [0, T]
nstep = [L / 0.02, T / 0.06]
[x,u] = heattheta(xspan, tspan, nstep, theta, alpha, u0, g, f);

PuntoMedido = 0.96;

position = find(x == PuntoMedido)
u(position)

if (u(position) < u0(PuntoMedido))
  printf("La temperatura en en el punto: %f ha bajado", PuntoMedido);
elseif (u(position) > u0(PuntoMedido))
  printf("La temperatura en el punto: %f ha subido", PuntoMedido);
else
  printf("La temperatura en el punto: %f se ha mantenido", PuntoMedido);
end


%Apartado d
a = 0
b = L = 2;
ua = 0 ; ub = uL = 0;
h = 0.02;
N = ((b - a) / h) - 1; %Numero Nodos internos
mu = alpha ; eta = 0; sigma = 0;
f = @(x) 2 - x

uex = @(x) (1 / alpha) * (((1 / 6) .* x^2) - x^2 + ((3/4) * x))

[x1, u1] = bvp(a, b, N, mu, eta, sigma, f, ua, ub)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% ERROR






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Newton