%% Part one: Calculation of eta
omega = 19;

S0 = @(x,y) cos(24*(x .^ 2 + y .^ 2)) .* exp(-900 * (x .^ 2 + y.^2));

etaf = @(X) S0(X(1, :), X(2, :)) .* cos(omega * vecnorm(X));

% we substitute R^2 with a square, S0 is small very fast possibly
% could try monte carlo integration with a normal distribution for better
% results for lower n (as seen on wikipedia)
n = 10^7;
s = 3; % e^(-2700) is prettei small, we ignore stuff outside of that (ish)
R2 = rand(2,n) * s;
muR2 = s^2;
eta = MonteCarlogeneric(etaf, R2, n, muR2);

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
% g = B.un
% gnoise=g+max(abs(g))*randn(size(g))*noiselevel;
IcIntegral = @(alpha) B.s(end) / length(B.s) * (vc(B.x, B.y, alpha)' * g);

IcShort = @(x0tilde, y0tilde, atilde, alpha) atilde*eta*vc(x0tilde, y0tilde, alpha);


nalphas = 30;
alphas = linspace(0,2*pi, nalphas);


tau = 1e-7;
F = @(x0tilde, y0tilde, atilde, alphas) atilde .* eta .* vc(x0tilde, y0tilde, alphas)' - IcIntegral(alphas);

x0tilde = ys + 0*rand();
y0tilde = xs + 0*rand();
a0tilde = a + 0*rand();

% for some reason this doesnt always seem to converge?? wtf
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
hold off;

%% 5
vs = @(x, y, alpha) sin(omega*(x*cos(alpha) + y*sin(alpha)));
IsIntegral = @(alpha) B.s(end) / length(B.s) * (vs(B.x, B.y, alpha)' * g);

Icprim = @(alpha, h) (IcIntegral(alpha + h) - IcIntegral(alpha - h))/(2*h);
Icprim0 = Icprim(0, 1e-10);
Icprim90deg = Icprim(pi/2, 1e-10);

x0 = Icprim90deg / (omega*IsIntegral(pi/2));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f is function, D is the set, muD is measure of D
% assumes f is vectorised properly
function I = MonteCarlogeneric(f, D, n, muD)
    I = muD / n * sum(f(D));
end




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

