function [Bound,Sol] = hhsolver(omega,S,N)

% Rummets dimensioner

Lx1 = 0.25;
Lx2 = 0.6;
Lx3 = 1;
Ly1 = 0.4;
Ly2 = 0.6;
Ly3 = 0.75;

% Diskretisering

dxmax = 0.05;  % Grövsta tillåtna diskretisering

N = round(dxmax*N)*20-1;
M = round(N*Ly3/Lx3);

dx = Lx3/(N+1);
dy = Ly3/(M+1);

x = dx:dx:Lx3-dx;
y = dy:dy:Ly3-dy;
[X,Y]=meshgrid(x,y);

% Ställ upp linjärt ekvationssystem

Ax = -2*diag(ones(N,1)) + diag(ones(N-1,1),1) + diag(ones(N-1,1),-1);
Ax = sparse(Ax)/dx^2;
Ay = -2*diag(ones(M,1)) + diag(ones(M-1,1),1) + diag(ones(M-1,1),-1);
Ay = sparse(Ay)/dy^2;

f = S(X,Y);

Ix = speye(N,N);
Iy = speye(M,M);
II = speye(M*N,M*N);

L = kron(Ax,Iy) + kron(Ix,Ay) + omega^2*II;
F = reshape(f,M*N,1);

% Reducera till ekvationer för inre punkter

w  = find(X<Lx1-dx/2 | (X>Lx2+dx/2 & Y<Ly2-dy/2) | Y<Ly1-dy/2);
L2 = L(w,w);
U2 = zeros(M*N,1);

% lös ekvationssystemet

U2(w)=L2\F(w);
U=reshape(U2,M,N);

% Lägg till kända randvärden till x, y, u

xf = 0:dx:Lx3;
yf = 0:dy:Ly3;
[Xf,Yf]=meshgrid(xf,yf);
Uf = [zeros(M+2,1) [zeros(1,N); U; zeros(1,N)] zeros(M+2,1)];

% Beräkna du/dn på randen

% Definition av rand-delarna
% binfo: [xstart xend ystart yend normalx normaly reverse]

binfo = [0 Lx3 0 0 0 -1 0; ...
    Lx3 Lx3 0 Ly2 1 0 0; ...
    Lx2 Lx3 Ly2 Ly2 0 1 1; ...
    Lx2 Lx2 Ly1 Ly2 -1 0 1; ...
    Lx1 Lx2 Ly1 Ly1 0 1 1; ...
    Lx1 Lx1 Ly1 Ly3 1 0 0; ...
    0 Lx1 Ly3 Ly3 0 1 1; ...
    0 0 0 Ly3 -1 0 1];

for j=1:8
    ww = find(Xf>binfo(j,1)-dx/2 & Xf<binfo(j,2)+dx/2 & Yf>binfo(j,3)-dy/2 & Yf<binfo(j,4)+dy/2);
    xj = Xf(ww); yj=Yf(ww);
    if (binfo(j,6)==0)  % Normal derivative in x dir
        sg = binfo(j,5);
        unj = (-2*Uf(ww-sg*(M+2))+0.5*Uf(ww-2*sg*(M+2)))/dy;
        ss = yj;
    else
        sg = binfo(j,6);
        unj = (-2*Uf(ww-sg)+0.5*Uf(ww-2*sg))/dy;
        ss = xj;
    end
    if (binfo(j,7) == 1)
        xj = flip(xj);
        yj = flip(yj);
        unj = flip(unj);
    end
    if (ss(end)<ss(1))
        ss = flip(ss);
    end
    
    if (j>1)
        UU(end) = (UU(end)+unj(1))/2;
        XX=[XX; xj(2:end)];  YY=[YY; yj(2:end)];  UU=[UU; unj(2:end)];  SS = [SS; SS(end)+ss(2:end)-ss(1)];
    else
        XX=xj;  YY=yj;  UU=unj;  SS = ss;
    end
    us = (UU(1)+UU(end))/2;
    UU(1)=us; UU(end)=us;
end

% Sätt alla värden utanför Omega till Not a Number.

wc = find(Xf>Lx1+dx/2 & (Xf<Lx2-dx/2 | Yf>Ly2+dy/2) & Yf>Ly1+dy/2);
Uf(wc)=0/0;
Xf(wc)=0/0;
Yf(wc)=0/0;

Sol=struct('x',Xf,'y',Yf,'u',Uf);
Bound=struct('s',SS,'x',XX,'y',YY,'un',UU);
%keyboard
%ev=eigs(-L2+omega^2*II(w,w),5,'smallestreal'); sqrt(ev)

end

