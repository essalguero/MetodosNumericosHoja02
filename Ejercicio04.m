clear all
clc

f = inline('0*x', 'x');

a = 0; b = 1;
h = 0.1
N = ((b - a) / h) - 1; %Numero Nodos internos
V1 = 1; V2 = 25; V3 = 100;
mu = 1;
eta = V1;
sigma = 0;

ua = 0; ub = 1; % Contorno

[x1, u1] = bvp(a, b, N, mu, V1, sigma, f, ua, ub)
[x2, u2] = bvp(a, b, N, mu, V2, sigma, f, ua, ub)
[x3, u3] = bvp(a, b, N, mu, V3, sigma, f, ua, ub)

sol1 = inline('(exp(V1 * x) - 1) ./ (exp(V1) - 1)', 'x')
sol2 = inline('(exp(V2 * x) - 1) ./ (exp(V2) - 1)', 'x')
sol3 = inline('(exp(V3 * x) - 1) ./ (exp(V3) - 1)', 'x')


subplot(1, 3, 1)
plot(x1, u1, x1, sol1(x1))
legend('aprox', 'exacta');
title('V = 1');

subplot(1, 3, 2)
plot(x2, u2, x2, sol2(x2))
legend('aprox', 'exacta');
title('V = 25');

subplot(1, 3, 3)
plot(x3, u3, x3, sol3(x3))
legend('aprox', 'exacta');
title('V = 100');