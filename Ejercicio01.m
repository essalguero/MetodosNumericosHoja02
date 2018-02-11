close all
clc

sol = @(t) 2 -exp(-10*t);

f = @(t,y) 10*(2 - y);
tspan = [0 3];
h1 = 0.3;
h2 = 0.15;
h3 = 0.05;
y0 = 1;

 

tspan = [0, 3];

[t1, y1] = feuler(f, tspan, y0, (tspan(2) - tspan(1))/h1);
[t2, y2] = feuler(f, tspan, y0, (tspan(2) - tspan(1))/h2);
[t3, y3] = feuler(f, tspan, y0, (tspan(2) - tspan(1))/h3);


figure
plot(t1, y1, t1, sol(t1))
title('h = 0.3');
legend('aprox', 'exacta');
figure
plot(t2, y2, t2, sol(t2))
title('h = 0.15');
legend('aprox', 'exacta');
figure
plot(t3, y3, t3, sol(t3))
title('h = 0.05');
legend('aprox', 'exacta');


h4 = 0.025;

[t4, y4] = beuler(f, tspan, y0, (tspan(2) - tspan(1))/h4);
figure
plot(t4, y4, t4, sol(t4))
title('Euler Implicito vs. Analitica');
legend('aprox', 'exacta');
[t4, y4](rows(t4),:)


[t5, y5] = cranknic(f, tspan, y0, (tspan(2) - tspan(1))/h4);
figure
plot(t5, y5, t5, sol(t5))
title('Crank-Nicholson vs. Analitica');
legend('aprox', 'exacta');
[t5', y5](rows(t5'),:)

valorAnalitico = sol(0.3)
valorBeuler = y4(lookup(t4, 0.3))
valorCN = y5(lookup(t5, 0.3))

errorBeuler = abs(valorAnalitico - valorBeuler)
errorCN = abs(valorAnalitico - valorCN)

if errorBeuler < errorCN
  disp("Beuler mas precios")
else
  disp("CN mas preciso")
end
