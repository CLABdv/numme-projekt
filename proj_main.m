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


xs  = 0.45;
ys  = 0.2;
a = 1;
S = @(x, y) a * S0(x-xs, y-ys);

[B, Sol] = hhsolver(omega, S, 300);
fprintf("We have chosen xs=%f, ys=%f and a=%f. Our value for eta is %f\n", xs, ys, a, eta)
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
%% Gauss Newton


% the vc(alpha)' * g makes each row dot product
% 4
g = B.un;
noiselevel = 0.2;
g=g+max(abs(g))*randn(size(g))*noiselevel; % gnoise
vs = @(x, y, alpha) sin(omega*(x*cos(alpha) + y*sin(alpha)));
vc = @(x, y, alpha) cos(omega*(x*cos(alpha) + y*sin(alpha)));
IcIntegral = @(alpha) B.s(end) / length(B.s) * (vc(B.x, B.y, alpha)' * g);
IsIntegral = @(alpha) B.s(end) / length(B.s) * (vs(B.x, B.y, alpha)' * g);


Icprim = @(alpha, h) (IcIntegral(alpha + h) - IcIntegral(alpha - h))/(2*h);
Icprim0 = Icprim(0, 1e-10);
Icprim90deg = Icprim(pi/2, 1e-10);


nalphas = 30;
alphas = linspace(0,2*pi, nalphas);


tau = 1e-7;
Ftemp = @(atilde, x0tilde, y0tilde, alphas) atilde .* eta .* vc(x0tilde, y0tilde, alphas)' - IcIntegral(alphas);

x0tilde = Icprim90deg / (omega*IsIntegral(pi/2));
y0tilde = -Icprim0 / (omega * IsIntegral(0));
a0tilde = 1 / eta * sqrt(IcIntegral(0)^2 + IsIntegral(0)^2);

X0 = [a0tilde; x0tilde; y0tilde];
F = @(X) Ftemp(X(1), X(2), X(3), alphas);
J = @(X) JacIc(X(1), X(2), X(3), alphas, eta, omega);
X = gaussNewton(X0, F, J,tau);

fprintf("With noiselevel = %f our estimations of a = %f, x0 = %f, y0 = %f\n", noiselevel, X(1), X(2), X(3));

%% 5 / 6
T = zeros(5,3);
for i = 1:5
    filename = sprintf("source%d.mat", i);
    load(filename);
    g=B.un;

    % these need to be redefined since our omega changes
    etaf = @(x ,y) S0(x, y) .* cos(omega * x);
    n = 10^7;
    s = 3;
    R2 = rand(2,n) * s - 1.5;
    fs = etaf(R2(1,:), R2(2,:));
    eta = s*s/n * sum(fs);

    vs = @(x, y, alpha) sin(omega*(x*cos(alpha) + y*sin(alpha)));
    vc = @(x, y, alpha) cos(omega*(x*cos(alpha) + y*sin(alpha)));
    IcIntegral = @(alpha) B.s(end) / length(B.s) * (vc(B.x, B.y, alpha)' * g);
    IsIntegral = @(alpha) B.s(end) / length(B.s) * (vs(B.x, B.y, alpha)' * g);
    Icprim = @(alpha, h) (IcIntegral(alpha + h) - IcIntegral(alpha - h))/(2*h);
    Icprim0 = Icprim(0, 1e-10);
    Icprim90deg = Icprim(pi/2, 1e-10);
        
    x0tilde = Icprim90deg / (omega*IsIntegral(pi/2));
    y0tilde = -Icprim0 / (omega * IsIntegral(0));
    a0tilde = 1 / eta * sqrt(IcIntegral(0).^2 + IsIntegral(0).^2);
    
    Ftemp = @(atilde, x0tilde, y0tilde, alphas) atilde .* eta .* vc(x0tilde, y0tilde, alphas)' - IcIntegral(alphas);
    F = @(X) Ftemp(X(1), X(2), X(3), alphas);
    J = @(X) JacIc(X(1), X(2), X(3), alphas, eta, omega);

    X0 = [a0tilde; x0tilde; y0tilde];
    X = gaussNewton(X0, F, J,tau);
    T(i, :) = X;
end
tab=array2table(T,'VariableNames',{'atilde', 'x0tilde', 'y0tilde'},'RowNames',{'source1','source2','source3','source4', 'source5'});
disp(tab);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%