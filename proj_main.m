%% Part one: Calculation of eta
omega = 19;

S0 = @(x,y) cos(24*(x .^ 2 + y .^ 2)) .* exp(-900 * (x .^ 2 + y.^2));

etaf = @(x ,y) S0(x, y) .* cos(omega * x);

% we substitute R^2 with a square, S0 is small very fast possibly
% could try monte carlo integration with a normal distribution for better
% results for lower n (as seen on wikipedia)
n = 10^7;
s = 3;
R2 = rand(2,n) * s - 1.5;

%eta = MonteCarlogeneric(etaf, R2, n, muR2);
fs = etaf(R2(1,:), R2(2,:));
eta = s*s/n * sum(fs);

%eta = integral2(etaf, -3, 3, -3, 3)

xs  = 0.45;
ys  = 0.2;
a = 1;
S = @(x, y) a * S0(x-xs, y-ys);

[B, Sol] = hhsolver(omega, S, 300);

%% Some random plotting
figure(1)
contour(Sol.x,Sol.y,Sol.u,20)
axis equal
hold on
plot(B.x,B.y,'k-','LineWidth',2)
[c,hnd]=contour(Sol.x,Sol.y,S(Sol.x,Sol.y),10); %Sol.S,10);
set(hnd,'Color','k','LineWidth',1.5)
hold off
axis off
%%

vc = @(x, y, alpha) cos(omega*(x*cos(alpha) + y*sin(alpha)));

% the vc(alpha)' * g makes each row dot product
% 4
g = B.un;
noiselevel = 0.5;
g=g+max(abs(g))*randn(size(g))*noiselevel; % gnoise
IcIntegral = @(alpha) B.s(end) / length(B.s) * (vc(B.x, B.y, alpha)' * g);

IcShort = @(x0tilde, y0tilde, atilde, alpha) atilde*eta*vc(x0tilde, y0tilde, alpha);


nalphas = 30;
alphas = linspace(0,2*pi, nalphas);


tau = 1e-7;
F = @(x0tilde, y0tilde, atilde, alphas) atilde .* eta .* vc(x0tilde, y0tilde, alphas)' - IcIntegral(alphas);

x0tilde = xs + 0.1*rand();
y0tilde = ys + 0.1*rand();
a0tilde = a + 0.1*rand();

while true
    dx =  -JacIc(x0tilde,y0tilde,a0tilde,alphas, eta) \ F(x0tilde, y0tilde, a0tilde, alphas);
    a0tilde = a0tilde + dx(1);
    x0tilde = x0tilde + dx(2);
    y0tilde = y0tilde + dx(3);
    if all (abs(dx) - tau < 0)
        break
    end
end

v = [a0tilde, x0tilde, y0tilde];

%% testplots

figure(2)
plot(alphas, IcShort(xs, ys, a, alphas));
hold on;
plot(alphas, IcIntegral(alphas));
plot(alphas, IcShort(x0tilde,y0tilde,a0tilde,alphas));


%% 5
vs = @(x, y, alpha) sin(omega*(x*cos(alpha) + y*sin(alpha)));
IsIntegral = @(alpha) B.s(end) / length(B.s) * (vs(B.x, B.y, alpha)' * g);

Icprim = @(alpha, h) (IcIntegral(alpha + h) - IcIntegral(alpha - h))/(2*h);
Icprim0 = Icprim(0, 1e-10);
Icprim90deg = Icprim(pi/2, 1e-10);

x0tilde = Icprim90deg / (omega*IsIntegral(pi/2));
y0tilde = -Icprim0 / (omega * IsIntegral(0));
a0tilde = 1 / eta * sqrt(IcIntegral(0)^2 + IsIntegral(0)^2);

while true
    dx =  -JacIc(x0tilde,y0tilde,a0tilde,alphas, eta) \ F(x0tilde, y0tilde, a0tilde, alphas);
    a0tilde = a0tilde + dx(1);
    x0tilde = x0tilde + dx(2);
    y0tilde = y0tilde + dx(3);
    if all (abs(dx) - tau < 0)
        break
    end
end

v = [a0tilde, x0tilde, y0tilde]

plot(alphas, IcShort(x0tilde,y0tilde,a0tilde,alphas));
hold off;




%% Optimal TV-placering

ys = 0.6;
xs = linspace(0.6, 1);
A = zeros(size(xs));
for i = 1:length(xs)
    S = @(x, y) a * S0(x-xs(i), y-ys);
    [B, Sol] = hhsolver(omega, S, 100);
    w = find(Sol.x<=0.25 & Sol.y>=0.5);
    A(i) = max(abs(Sol.u(w)))/max(abs(Sol.u(:)));
end
figure(4);
plot(xs, A);
xs = 0.77; % approximative global minimum
%% more exact solution



tau = 1e-4;
eps = 0.02;

a = xs - eps;
b = xs + eps;

y = (sqrt(5) - 1)/2 * (b-a) + a;
x = b - (sqrt(5) - 1)/2 * (b-a);

f1 = f(x, ys, 1, S0, 1000);
f2 = f(y, ys, 1, S0, 1000);

while b - a > tau
    if f1 < f2
        a = x;
        x = b - (sqrt(5) - 1)/2 * (b-a);
        f1 = f(x, ys, 1, S0, 1000);
    else 
        b = y;
        y = (sqrt(5) - 1)/2 * (b-a) + a;
        f2 = f(y, ys, 1, S0, 1000);
    end
    disp(b-a)
end

disp("finshid")

%% Grafer av bäst placering
xs = a;
S = @(x, y) 1 * S0(x-xs, y-ys);
[B, Sol] = hhsolver(omega, S, 1000);
figure(5)
mesh(Sol.x,Sol.y,Sol.u)

figure(6)
contour(Sol.x,Sol.y,Sol.u,20)
axis equal
hold on
plot(B.x,B.y,'k-','LineWidth',2)
[c,hnd]=contour(Sol.x,Sol.y,S(Sol.x,Sol.y),10); %Sol.S,10);
set(hnd,'Color','k','LineWidth',1.5)
hold off
axis off

%% Grafer av relativt dålig placering (please consult the graph)
xs = 0.87;
S = @(x, y) 1 * S0(x-xs, y-ys);
[B, Sol] = hhsolver(omega, S, 1000);
figure(7)
mesh(Sol.x,Sol.y,Sol.u)

figure(8)
contour(Sol.x,Sol.y,Sol.u,20)
axis equal
hold on
plot(B.x,B.y,'k-','LineWidth',2)
[c,hnd]=contour(Sol.x,Sol.y,S(Sol.x,Sol.y),10); %Sol.S,10);
set(hnd,'Color','k','LineWidth',1.5)
hold off
axis off

%% Two dimensional optimisation
% Grovlokalisering med n = 100
npoints = 20;

xs = linspace(0.6, 1, npoints);
ys = linspace(0, 1, npoints);

As = zeros(npoints, npoints);

for i = 1:length(xs)
    for j = 1:length(ys)
        As(i,j) = f(xs(i), ys(j), 1, S0, 100);
    end
    disp(i)
end

figure(9)
mesh(xs, ys, As)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Jacobi matrix for Ic
function J = JacIc(x0tilde, y0tilde, atilde, alphas, eta)
     omega = 19;
     J = zeros(length(alphas),3);
     for i = 1:length(alphas)
         J(i,1) = eta * cos(omega * (x0tilde .* cos(alphas(i)) + y0tilde .* sin(alphas(i)))); %dF/da
         J(i,2) = -atilde * eta .* sin(omega *(x0tilde .*cos(alphas(i)) + y0tilde .*sin(alphas(i))))*omega .*cos(alphas(i)); %dF/dx0
         J(i,3) = -atilde * eta .* sin(omega *(x0tilde .*cos(alphas(i)) + y0tilde .*sin(alphas(i))))*omega .*sin(alphas(i)); %dF/dy0
     end
end


function A = f(xs, ys, a, S0, n)
    omega = 19;
    S = @(x, y) a * S0(x-xs, y-ys);
    [~, Sol] = hhsolver(omega, S, n);
    w = find(Sol.x<=0.25 & Sol.y>=0.5);
    A = max(abs(Sol.u(w)))/max(abs(Sol.u(:)));
end