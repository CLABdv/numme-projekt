%% Part one: Calculation of eta
omega = 19;

S0 = @(x,y) cos(24*(x .^ 2 + y .^ 2)) .* exp(-900 * (x .^ 2 + y.^2));
etafinternal = @(X) S0(X(1, :), X(2, :)) .* cos(omega * vecnorm(X));
etaf = @(X) etafinternal(X);

% we substitute R^2 with a square, S0 is small very fast possibly
% could try monte carlo integration with a normal distribution for better
% results for lower n (as seen on wikipedia)
n = 10^7;
s = 3; % e^(-2700) is prettei small, we ignore stuff outside of that (ish)
R2 = rand(2,n) * s;
muR2 = s^2;
eta = MonteCarlogeneric(etaf, R2, size(R2, 2), muR2);

xs  = 0.45;
ys  = 0.2;
S = @(x, y) S0(x-xs, y-ys);
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

g = B.un;
vc = @(x, y, alpha) cos(omega*(cos(alpha) .* x + sin(alpha) .* y));


% the vc(alpha)' * g should be summing all the points in the integrand
IcIntegral = @(alpha) B.s(end) / size(B.s,1) * (vc(B.x, B.y, alpha)' * g);
IcShort = @(x0tilde, y0tilde, atilde, alpha) a*eta*vc(x0tilde, y0tilde, alpha);

%alphas = linspace(0,2*pi, nalphas);
alpha = pi;
Y = @(x0tilde, y0tilde, atilde, alphas) atilde * eta * vc(x0tilde, y0tilde, alphas) - IcIntegral(alphas);

x0 = xs + 0.5*rand();
y0 = ys + 0.5*rand();
a = 2;
tau = 1e-10;
while true
    dx = -Y(x0, y0, a, alpha) \ JacIc(x0,y0,a,alpha, eta);
    if all (abs(dx) - tau < 0)
        break
    end
    x0 = x0 + dx(1)
    y0 = y0 + dx(2)
    a = a + dx(3)
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f is function, D is the set, muD is measure of D
% assumes f is vectorised properly
function I = MonteCarlogeneric(f, D, n, muD)
    I = muD / n * sum(f(D));
end

% Jacobi matrix for Ic (well a gradient cuz real valued but y'know)
function J = JacIc(x0tilde, y0tilde, atilde, alpha, eta)
     omega = 19;
     J = zeros(1,3);
     J(1) = eta * cos(omega * (x0tilde .* cos(alpha) + y0tilde .* sin(alpha)));
     J(2) = -atilde * eta .* sin(omega *(x0tilde .*cos(alpha) + y0tilde .*sin(alpha)))*omega .*cos(alpha);
     J(3) = -atilde * eta .* sin(omega *(x0tilde .*cos(alpha) + y0tilde .*sin(alpha)))*omega .*sin(alpha);
end

