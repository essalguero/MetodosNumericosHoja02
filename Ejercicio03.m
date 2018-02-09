clear all
clc




S = 100;
D = 8.8 * 10^7;
q = 200;
l = 50;

a = 0;
b = l;

Wa = 0;
Wl = 0;

mu = 1;
eta = 0;
sigma = -S / D;

%f = inline('-(((q * l)/(2 * D)) * x) + (((q)/(2 * D)) * x .^2)', 'x')

f = @(x) ((-q * l) / (2 * D)) .* x + ((q / (2 * D)) .* x.^2);


h = 1;

N = ((l - 0) / h) - 1 %Numero Nodos internos

[x1, u1] = bvp(a, b, N, mu, eta, sigma, f, Wa, Wl)

% Hay partes de la placa en la que la deflexion es mayor de abs(1)

