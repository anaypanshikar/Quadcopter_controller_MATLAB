syms s kp1 tau
kp1 = 46.34;
G = tf([kp1],[1,0,0]);
[A1,B1,C1,D1] = tf2ss([kp1],[1,0,0]);
sys = ss(A1,B1,C1,D1);
c = canon(sys,'controller');
global A B C D Am Bm Cm Dm
A = c.A;
B = c.B;
C = c.C;
D = c.D;

%Model poles
z = 0.3;
Ts = 1.2;
wn = 4/(z*Ts);
wd = wn*sqrt(1-z*z);
sigma = 4/Ts;
r1 = -sigma+j*wd;
r2 = -sigma-j*wd;
num1 = [wn^2];
den1 = sym2poly(expand((s-r1)*(s-r2)));
Gm = tf(num1,den1);
%step(Gm)

m = canon(Gm,'controller');
Am = m.A;
Bm = m.B;
Cm = m.C;
Dm = m.D;

x0 = [0;0];
x0m = [0;0];

P = sylvester(Am',Am,-1*eye(2));
%eig(P)

syms uc
Lr = [0.15];
L = [1.9,9.1];
%vpa(A-B*L,4);
%vpa(B*Lr,4);

x = x0; %initial state
xm = x0m;
uc = 1; %reference input signal

times = [0:0.02:5];
figure
hold on
for t = times
    u = sysinput(Lr,L,x,uc);
    [x,xm] = update(x,xm,u,uc);
    e = error(x,xm);
    [L,Lr] = paraupdate(L,Lr,P,e,x,uc);
    scatter(t,C*x,15,'filled','red')
    hold on
    scatter(t,Cm*xm,15,'filled','blue')
    ylim([-2,3])
    xlim([0,5])
    legend('System', 'Model')
    pause(0.00001)
end

function [x,xm] = update(x,xm,u,uc)
    global A B Am Bm
    dt = 0.02;
    xdot = A*x + B*u;
    xmdot = Am*xm + Bm*uc;
    x = x + xdot*dt;
    xm = xm + xmdot*dt;
end

function u = sysinput(Lr,L,x,uc)
    u = Lr*uc - L*x;
end

function [L,Lr] = paraupdate(L,Lr,P,e,x,uc)
    dt = 0.02;
    temp1 = P*e*x';
    c123dot = temp1(1,:);
    temp2 = -1*P*e*uc;
    c4dot = temp2(1);
    L = L + c123dot*dt;
    Lr = Lr + c4dot*dt;
end

function e = error(x,xm)
    e = x-xm;
end
