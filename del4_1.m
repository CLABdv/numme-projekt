%% Optimal TV-placering
disp("Optimal TV-placering:")
omega = 30;
fprintf("omega = %f\n", omega)

% kod för att göra grovsökning efter minimum längs väggen
ys = 0.6;
xs = linspace(0.6, 1);
A = zeros(size(xs));
for i = 1:length(xs)
    S = @(x, y) a * S0(x-xs(i), y-ys);
    [B, Sol] = hhsolver(omega, S, 100);
    w = find(Sol.x<=0.25 & Sol.y>=0.5);
    A(i) = max(abs(Sol.u(w)))/max(abs(Sol.u(:)));
end
figure(3);
plot(xs, A);
legend("A", 'Location', 'southeast')
xs = 0.67; % approximative global minimum
%% grejer för att ge graf för lösningen

xs = a;
S = @(x, y) 1 * S0(x-xs, y-ys);
[B, Sol] = hhsolver(omega, S, 1000);