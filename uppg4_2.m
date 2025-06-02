%% Grafer av relativt dålig placering
xs = 0.75;

S = @(x, y) 1 * S0(x-xs, y-ys);
[B, Sol] = hhsolver(omega, S, 1000);