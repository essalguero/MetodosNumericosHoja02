clear all
clc

format long

sol = @(t) exp(((1/2) * t .^ 2) - t) ./ (1 + exp(((1/2) * t .^ 2) - t))

f = @(t, y) (t - 1).*(y - y.^2)
y0 = 0.5
tspan = [0 2];
h1 = 0.004;

[t1, y1] = rk2(f, tspan, y0, (tspan(2) - tspan(1))/h1);
[t2, y2] = rk3(f, tspan, y0, (tspan(2) - tspan(1))/h1);

y1(rows(t1))
y2(rows(t2))

valorAnalitico = sol(1.24)
valorRk2 = y1(lookup(t1, 1.24))
valorRk3 = y2(lookup(t2, 1.24))

errorRk2 = abs(valorAnalitico - valorRk2)
errorRk3 = abs(valorAnalitico - valorRk3)

if errorRk2 < errorRk3
  disp("errorRk2 mas precios")
else
  disp("errorRk3 mas preciso")
end


%Calcular las soluciones exactas
soluciones = zeros(rows(t1), 2);
for i = 1:rows(t1)
  soluciones(i,:) = [t1(i), sol(t1(i), 1)];
end

soluciones(rows(t1))

erroresRk2 = abs(y1 .- soluciones(:,2));

erroresRk3 = abs(y2 .- soluciones(:,2));

figure
plot(t1, erroresRk2)
title("Error")
legend("Error errorRk2")

figure
plot(t2, erroresRk3)
title("Error")
legend("Error errorRk3")

tablaTiempos04 = find(soluciones(:,2) < 0.4);
soluciones(min(tablaTiempos04), 1)
soluciones(max(tablaTiempos04), 1)


I = [0:10^-6:2];
soluciones = sol(I);
minimo = min(find(soluciones < 0.4)) * 10^-6
maximo = max(find(soluciones < 0.4)) * 10^-6

figure
plot(I, sol(I))
title("Funcion original")
legend("Sol")

minimoRK2 = min(find(y1 < 0.4)) * h1
maximoRK2 = max(find(y1 < 0.4)) * h1

figure
plot(I, sol(I))
title("Funcion RK2")
legend("RK2")

minimoRK3 = min(find(y2 < 0.4)) * h1
maximoRK3 = max(find(y2 < 0.4)) * h1

figure
plot(I, sol(I))
title("Funcion RK3")
legend("RK3")



% bisection