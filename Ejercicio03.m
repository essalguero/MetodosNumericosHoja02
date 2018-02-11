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

f = @(x) ((-q * l) / (2 * D)) .* x + ((q / (2 * D)) .* x.^2);


h = 1;

N = ((l - 0) / h) - 1 %Numero Nodos internos

[x1, u1] = bvp(a, b, N, mu, eta, sigma, f, Wa, Wl);

fprintf("Hay partes de la placa en la que la abs(deflexion) es mayor de 0.1\n\n");

printf("posicion - valor\n");
Posiciones = find(abs(u1) > 0.1);

[Posiciones, u1(Posiciones)]