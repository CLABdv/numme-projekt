%% Part one: Calculation of eta

omega = 19;

S0 = @(x,y) cos(24*(x .^ 2 + y .^ 2)) .* exp(-900 * (x .^ 2 + y.^2));
etafinternal = @(X) S0(X(1, :), X(2, :)) .* cos(omega * vecnorm(X));
etaf = @(X) etafinternal(X);

% we substitute R^2 with a square, S0 is small very fast possibly
% could try monte carlo integration with a normal distribution for better
% results for lower n (as seen on wikipedia)
n = 10^7;
s = 3; % e^(-2700) is prettei small
R2 = rand(2,n) * s;
muR2 = s^2;
eta = MonteCarlogeneric(etaf, R2, size(R2, 2), muR2);

xs  = 0.45;
ys  = 0.2;
S = @(x, y) S0(x-xs, y-ys);
[Bound, Sol] = hhsolver(omega, S, 300);

%%
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

g = Bound.un






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f is function, D is the set, muD is measure of D
% assumes f is vectorised properly
function I = MonteCarlogeneric(f, D, n, muD)
    I = muD / n * sum(f(D));
end