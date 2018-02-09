clear all
clc


T = 30; L = 2; alpha = 0.004

theta = 0; %Euler Explicito

f = @(x, t) 2 - x + exp(-2 * t)
g = @(x, t) (2 / (1 + (t - (1 / 2))^2)) - x.*(2 - x)

u0 = inline('2.*(1 - x).^2.*(heaviside(x-1/2)-heaviside(x-3/2))', 'x')

xspan = [0, L]
tspan = [0, T]
nstep = [L / 0.02, T / 0.06]
heattheta(xspan, tspan, nstep, theta, alpha, u0, g, f)