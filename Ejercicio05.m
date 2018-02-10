clear
clc

format long

f = @(x, y) 8 * pi^2 * sin(2 * pi .* x) * cos(2 * pi .* y)
DeltaU = f;
h1 = 1 / 10;
h2 = 1 / 20;

a = 0; b = 1; c = 0; d = 1;

bound = @(x) sin(2 * pi * x)

uex = @(x, y) sin(2 * pi * x) * cos(2 * pi * y)

[u1, x1, y1, error1] = poissonfd(a, c, b, d, (b - a) / h1, (b - a) / h1, f, bound, uex);
[u2, x2, y2, error2] = poissonfd(a, c, b, d, (b - a) / h2, (b - a) / h2, f, bound, uex);

u1
u2

printf("Error1: %f\n", error1);
printf("Error2: %f\n", error2);

mesh(x1, y1, u1);
mesh(x2, y2, u2);

uex(x1, y1) - error1