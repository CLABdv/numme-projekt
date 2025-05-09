%% Two dimensional optimisation
% Grovlokalisering med n = 100
npoints = 20;

xs = linspace(0.6, 1, npoints);
ys = linspace(0, 1, npoints);

As = zeros(npoints, npoints);

xsMin = NaN;
ysMin = NaN;
AMin = inf;

for i = 1:length(xs)
    for j = 1:length(ys)
        As(i,j) = f(xs(i), ys(j), 1, S0, 100);
        if As(i, j) < AMin
            xsMin = xs(i);
            ysMin=ys(j);
            AMin = As(i,j);
        end
    end
end

figure(8)
mesh(xs, ys, As)
disp("grov minimum [x;y]:")
disp([xsMin; ysMin])
g = @(X) f(X(1), X(2), 1, S0, 1000);
options = optimset('Display', 'iter');
x = fminsearch(g, [xsMin;ysMin], options);

S = @(x, y) 1 * S0(x-xsMin, y-ysMin);
[B, Sol] = hhsolver(omega, S, 1000);