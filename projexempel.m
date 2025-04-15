omega = 30;
xs  = 0.45;
ys  = 0.2;
aa  = 1;

%omega = 19;
%xs  = 0.8;
%ys  = 0.43;
%aa  = 1;

S0 = @(x,y) exp(-400*(x.^2+y.^2));

S = @(x,y) aa*S0(x-xs,y-ys);

[B,Sol]=hhsolver(omega,S,200); 

figure(1)

surf(Sol.x,Sol.y,Sol.u)
view(0,90)
shading interp
axis equal

figure(2)

surf(Sol.x,Sol.y,Sol.u)
shading interp
light                    
lighting gouraud       % Tar lång tid när N är stor
view(-60,60)

figure(3)
%surf(Sol.x,Sol.y,Sol.S)
surf(Sol.x,Sol.y,S(Sol.x,Sol.y))
shading interp
light
lighting gouraud       % Tar lång tid när N är stor
colormap('autumn')
material shiny
view(-60,60)

figure(4)
contour(Sol.x,Sol.y,Sol.u,20)
axis equal
hold on
plot(B.x,B.y,'k-','LineWidth',2)
[c,hnd]=contour(Sol.x,Sol.y,S(Sol.x,Sol.y),10); %Sol.S,10);
set(hnd,'Color','k','LineWidth',1.5)
hold off
axis off

figure(5)

plot3(B.x,B.y,B.un,'b-','LineWidth',2)
hold on
contour(Sol.x,Sol.y,Sol.u,20)
plot3(B.x,B.y,zeros(size(B.x)),'k-','LineWidth',2)
hold off

figure(6)
mesh(Sol.x,Sol.y,Sol.u)



