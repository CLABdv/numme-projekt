proj_main;
%% Placering av TV längs väggen
del4_1;

xs = 0.67; % startgissning utifrån grovsearch

tau = 1e-4;
eps = 0.02;

a = xs - eps;
b = xs + eps;

func = @(x) f(x, ys, 1, S0, 1000);

[lower, upper] = phi_src(func, a, b, tau);


fprintf("TV'n har sin optimala placering x0 \\in [%f;%f] när y är fixerad på väggen, den relativa ljudnivån blir då %f\n", lower, upper, f(lower, ys, 1, S0, 1000));

% Coola grafer för att se hur ljudet blir
figure(4)
mesh(Sol.x,Sol.y,Sol.u)

figure(5)
contour(Sol.x,Sol.y,Sol.u,20)
axis equal
hold on
plot(B.x,B.y,'k-','LineWidth',2)
[c,hnd]=contour(Sol.x,Sol.y,S(Sol.x,Sol.y),10); %Sol.S,10);
set(hnd,'Color','k','LineWidth',1.5)
hold off
axis off
%% Sämre placering för jämförelse
uppg4_2;

fprintf("x0=%f är en relativt dålig placering, relativ ljudnivå A = %f\n", xs, f(xs, ys, 1, S0, 1000))

figure(6)
mesh(Sol.x,Sol.y,Sol.u)

figure(7)
contour(Sol.x,Sol.y,Sol.u,20)
axis equal
hold on
plot(B.x,B.y,'k-','LineWidth',2)
[c,hnd]=contour(Sol.x,Sol.y,S(Sol.x,Sol.y),10); %Sol.S,10);
set(hnd,'Color','k','LineWidth',1.5)
hold off
axis off
%% Tvådimensionella fallet
del4_3;

fprintf("Den optimala placeringen av TV'n i det högra rummet är [%f, %f], då är den relativa ljudnivån A = %f\n", x(1), x(2), f(x(1),x(2),1,S0,1000))

% mer grafer

figure(9)
mesh(Sol.x,Sol.y,Sol.u)

figure(10)
contour(Sol.x,Sol.y,Sol.u,20)
axis equal
hold on
plot(B.x,B.y,'k-','LineWidth',2)
[c,hnd]=contour(Sol.x,Sol.y,S(Sol.x,Sol.y),10); %Sol.S,10);
set(hnd,'Color','k','LineWidth',1.5)
hold off
axis off

%%