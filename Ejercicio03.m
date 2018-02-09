clear all
clc




S = 100;
D = 8.8 * 10^7;
q = 200;
l = 50;

WInit = 0;
WFin = 1;

mu = 1;
eta = 0;
sigma = -S / D;

f = inline('-(((q * l)/(2 * D)) * x) + (((q)/(2 * D)) * x ^2)', 'x')

h = 1

N = ((WFin - WInit) / h) - 1 %Numero Nodos internos

[x1, u1] = bvp(WInit, WFin, N, mu, eta, sigma, f, 0, 1)