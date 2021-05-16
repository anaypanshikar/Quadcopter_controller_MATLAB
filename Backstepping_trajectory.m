close all
clear all 
clc

global M g Ixx Iyy Izz l b d
M = 0.8;
g = 9.81;
l = 0.25;
b = 3*10^-6;
d = 10^-7;
Ixx = 0.005;
Iyy = 0.005;
Izz = 0.01;

I = [Ixx 0 0; 0 Iyy 0; 0 0 Izz];

%% initial and desired parameters
euler = [-10*pi/180 5*pi/180, 10*pi/180]';
eulerdot = [0 0 0]';
pos = [1 0 0]';
posdot = [0 0 0]';

psid = -10*pi/180;
% xd = 5;
% yd = 4;
% zd = 8;

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
c8 = 31.5;
c9 = 1.725;
c10 = 0.63;
c11 = 1.725;
c12 = 0.63;

u1 = M*g;
u2 = 0;
u3 = 0;
u4 = 0;
ux = 0;
uy = 0;

dt = 0.01;
for k = 1:2000
    
    xdes = [0;0;0;0;psid;0;0.1*k*dt;0.1;cos(k*dt);-sin(k*dt);sin(k*dt);cos(k*dt)];
    xddes = [0;0;0;0;0;0;0;0;-sin(k*dt);-cos(k*dt);cos(k*dt);-sin(k*dt)];
    xdddes = [0;0;0;0;0;0;0;0;-cos(k*dt);sin(k*dt);-sin(k*dt);-cos(k*dt)];
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
figure(1)
plot(t,pos(1,:),t,pos(2,:),t,pos(3,:));
legend("x","y","z")
xlabel("time (s)")
ylabel("position (m)")

%numel(t)
figure(2)
plot(t,euler(1,:),t,euler(2,:),t,euler(3,:));
legend("phi","theta","psi")
xlabel("time (s)")
ylabel("euler angles (degrees)")
% 
% %numel(t)
% figure(3)
% plot(t,eulerdot(1,:)*180/pi,t,eulerdot(2,:)*180/pi,t,eulerdot(3,:)*180/pi);
% legend("phidot","thetadot","psidot")
% xlabel("time (s)")
% ylabel("anglular velocities (degrees/s)")
% 
t = linspace(0,20,2001);
%numel(t)
figure(4)
plot(t,u2*-0.25,t,u3*-0.25,t,u4*-1);
legend("u2","u3","u4")
xlabel("time (s)")
ylabel("torque control input")
%ylim([-0.05 0.15])

t = linspace(0,20,2001);
%numel(t)
figure(5)
plot(t,u1);
legend("u1")
xlabel("time (s)")
ylabel("thrust control input")
%ylim([4.5 6])

figure(6)
plot3(pos(1,:),pos(2,:),pos(3,:))
hold on
plot3(ref(1,:),ref(2,:),ref(3,:))
xlabel("X")
ylabel("Y")
zlabel("Z")
zlim([0 2.5])
hold off
grid on
legend('Actual Path','Reference Path')

t = linspace(0,20,2000);
figure(7)
plot(t,error_x(1,:),t,error_x(2,:),t,error_x(3,:)),legend('ex','ey','ez')
xlabel("time")
ylabel("error (m)")
