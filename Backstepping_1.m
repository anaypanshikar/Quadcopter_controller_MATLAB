close all
clear all 
clc

global M g Ixx Iyy Izz l b d
M = 0.5;
g = 9.81;
l = 0.25;
b = 3*10^-6;
d = 10^-7;
Ixx = 0.005;
Iyy = 0.005;
Izz = 0.01;

I = [Ixx 0 0; 0 Iyy 0; 0 0 Izz];

%% initial and desired parameters
euler = [10*pi/180 -10*pi/180, 0*pi/180]';
eulerdot = [0 0 0]';
pos = [0 0 1]';
posdot = [0 0 0]';

psid = 10*pi/180;
xd = 5;
yd = 6;
zd = 4;

%% Simulation

a1 = (Iyy-Izz)/Ixx;
a2 = (Izz-Ixx)/Iyy;
a3 = (Ixx-Iyy)/Izz;
b1 = l/Ixx;
b2 = l/Iyy;
b3 = 1/Izz;

x = [euler(1);eulerdot(1);euler(2);eulerdot(2);euler(3);eulerdot(3);pos(3);posdot(3);pos(1);posdot(1);pos(2);posdot(2)];

c1 = 6;
c2 = 4;
c3 = 5.2;
c4 = 3;
c5 = 7;
c6 = 4.9;
c7 = 6.125;
c8 = 10.5;
c9 = 0.1225;
c10 = 1.63;
c11 = 0.1225;
c12 = 1.63;

u1 = M*g;
u2 = 0;
u3 = 0;
u4 = 0;
ux = 0;
uy = 0;

dt = 0.01;
for k = 1:2000
    
    xdes = [0;0;0;0;psid;0;zd;0;xd;0;yd;0];
    xddes = [0;0;0;0;0;0;0;0;0;0;0;0];
    xdddes = [0;0;0;0;0;0;0;0;0;0;0;0];
    ref(:,k) = [xdes(9),xdes(11),xdes(7)]';
    
    e1 = xdes(1)-x(1,k);
    e2 = x(2,k)-xddes(1)-c1*e1;
    e3 = xdes(3)-x(3,k);
    e4 = x(4,k)-xddes(3)-c3*e3;
    e5 = xdes(5)-x(5,k);
    e6 = x(6,k)-xddes(5)-c5*e5;
    e7 = xdes(7)-x(7,k);
    e8 = x(8,k)-xddes(7)-c7*e7;
    e9 = xdes(9)-x(9,k);
    e10 = x(10,k)-xddes(9)-c9*e9;
    e11 = xdes(11)-x(11,k);
    e12 = x(12,k)-xddes(11)-c11*e11;
    
    u2(k+1) = (e1 - a1*x(4,k)*x(6,k) + xdddes(1) - c2*e2 + c1*(xddes(1)-x(2,k)))./b1;
    u3(k+1) = (e3 - a2*x(2,k)*x(6,k) + xdddes(3) - c4*e4 + c3*(xddes(3)-x(4,k)))./b2;
    u4(k+1) = (e5 - a3*x(2,k)*x(4,k) + xdddes(5) - c6*e6 + c5*(xddes(5)-x(6,k)))./b3;
    u1(k+1) = M*(g + e7 + xdddes(7) - c8*e8 + c7*(xddes(7)-x(8,k)))./(cos(x(1,k))*cos(x(3,k)));
    
    ux(k+1) = M*(e9 + xdddes(9) - c10*e10 + c9*(xddes(9)-x(10,k)))./u1(k);
    uy(k+1) = M*(e11 + xdddes(11) - c12*e12 + c11*(xddes(11)-x(12,k)))./u1(k);
    
    xdot = [x(2,k); x(4,k)*x(6,k)*a1 + u2(k)*b1; x(4,k); x(2,k)*x(6,k)*a2 + u3(k)*b2; x(6,k); x(2,k)*x(4,k)*a3 + u4(k)*b3; x(8,k); u1(k)/M*cos(x(1,k))*cos(x(3,k)) - g; x(10,k); u1(k)*ux(k)/M; x(12,k); u1(k)*uy(k)/M];
    
    x(:,k+1) = x(:,k) + xdot*dt;
    
    euler(:,k+1) = [x(1,k+1);x(3,k+1);x(5,k+1)];
    eulerdot(:,k+1) = [x(2,k+1);x(4,k+1);x(6,k+1)];
    pos(:,k+1) = [x(9,k+1);x(11,k+1);x(7,k+1)];
    posdot(:,k+1) = [x(10,k+1);x(12,k+1);x(8,k+1)];
    
    error_x(:,k) = ref(:,k) - pos(:,k);
end

%% Plots

t = 0:0.01:20;
%numel(t)
figure
plot(t,pos(1,:),t,pos(2,:),t,pos(3,:));
legend("x","y","z")
xlabel("time (s)")
ylabel("position (m)")
ylim([0 8])

t = 0:0.01:10;
%numel(t)
figure
plot(t,euler(1,1:1001)*180/pi,t,euler(2,1:1001)*180/pi,t,euler(3,1:1001)*180/pi);
legend("\phi","\theta","\psi")
xlabel("time (s)")
ylabel("euler angles (degrees)")
ylim([-15 15])

t = 0:0.01:10;
%numel(t)
figure
plot(t,eulerdot(1,1:1001)*180/pi,t,eulerdot(2,1:1001)*180/pi,t,eulerdot(3,1:1001)*180/pi);
legend("$\dot{\phi}$","$\dot{\theta}$","$\dot{\psi}$",'interpreter','latex')
xlabel("time (s)")
ylabel("angular velocities (deg/s)")

t = linspace(0,20,2001);
%numel(t)
figure
plot(t,u2*-0.25,t,u3*-0.25,t,u4*-1);
legend("U2","U3","U4")
xlabel("time (s)")
ylabel("Torques (Nm)")
%ylim([-0.05 0.15])

t = linspace(0,20,2001);
%numel(t)
figure
plot(t,u1);
legend("u1")
xlabel("time (s)")
ylabel("thrust control input")
%ylim([4.5 6])

t = linspace(0,20,2000);
figure
plot(t,error_x(1,:),t,error_x(2,:),t,error_x(3,:)),legend('ex',"ey","ez")
xlabel('time')
ylabel('Error (m)')
