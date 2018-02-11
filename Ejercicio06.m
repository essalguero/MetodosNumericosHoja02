clear all
clc

format long

% Apartado a
printf("\n\n\nApartado a\n\n");

T = 30; L = 2; alpha = 0.004

theta = 1/2; %Euler Explicito

f = @(x, t) 2 .- x + exp(-2 .* t)
g = @(x, t) (2 / (1 .+ ((t .- 1) / 2).^2)) - x.*(2 .- x)

u0 = inline('2.*(1 - x).^2.*(heaviside(x-1/2)-heaviside(x-3/2))', 'x')

xspan = [0, L];
tspan = [0, T];
nstep = [L / 0.02, T / 0.06];
[x, u] = heattheta(xspan, tspan, nstep, theta, alpha, u0, g, f);

% Apartado b
printf("\n\n\nApartado b\n\n");
theta = 1;
nstep = [L / (0.02 * 2), T / (0.06 * 2)];
[xEE, uEE] = heattheta(xspan, tspan, nstep, theta, alpha, u0, g, f);
MaxTemp = max(uEE)
MinTemp = min(uEE)

% Apartado c

printf("\n\n\nApartado c\n\n");
PuntoMedido = 0.96;
TemperaturaInicial = u0(PuntoMedido)

posicionX = lookup(x, PuntoMedido);
TemperaturaFinal = u(posicionX)

printf('Al principio, la temperatura en %f era de %f\n', PuntoMedido, TemperaturaInicial);
printf('Despues es de %f\n', TemperaturaFinal);
printf('Se ha calentado\n\n');


%Apartado d

printf("\n\n\nApartado d\n\n");
printf("\n\n\nApartado d.i\n\n");

a = 0
b = L = 2;
ua = 0 ; ub = uL = 0;
h = 0.02;
N = ((b - a) / h) - 1; %Numero Nodos internos
mu = alpha ; eta = 0; sigma = 0;
f = @(x) 2 - x

uex = @(x) (1 / alpha) * (((1 / 6) .* x.^3) - x.^2 + ((4/3) .* x))

[x1, u1] = bvp(a, b, N, mu, eta, sigma, f, ua, ub);

MaxTemp = max(u1)
PositionMaxTemp = find(u1 > 128.296)
XValue = x1(PositionMaxTemp)
ValorExacto = uex(XValue)

ErrorCometido = abs(ValorExacto - MaxTemp)



printf("\n\n\nApartado d.ii\n\n");
% Newton
tol = 10^-3;
maxIter = 1000;
x0_75 = 0.3;
x1_75 = 1.5;

uex75 = @(x) (((1 / alpha) * (((1 / 6) .* x.^3) - x.^2 + ((4/3) .* x))) - 75)

df = @(x) (1/alpha) * (((1/2) .* x.^2) - (2 .* x) + (4/3));

[zero1, res1, niter1] = newton(uex75, df, x0_75, tol, maxIter);
[zero2, res2, niter2] = newton(uex75, df, x1_75, tol, maxIter);

rangoNewton = [zero1, zero2]


printf("\n\n\nApartado d.iii\n\n");

phi = @(x) x - (uex75(x) / df(x));

[xPuntoFijo, niter] = puntofijo2(phi, rangoNewton, tol, maxIter)

printf("Diferencia entre Newton y PuntoFijo\n");
xPuntoFijo - rangoNewton