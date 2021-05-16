close all
clear all 
clc 
% physical constants     
g = 9.81; 
M = 0.8; 
L = 0.25;
Ixx = 5*10^-3;
Iyy = 5*10^-3;
Izz = 10*10^-3; 
k = 3*10^-6; 
b = 10^-7;
I = diag([Ixx Iyy Izz]);

%constants
Kv = [0.63 ; 0.63; 31.5];  %for positions
Kp = [1.725 ; 1.725; 6.125];
Kpe = [6; 5.2; 7];      %for eular angles
Kpw = [4; 3; 4.9];       % for w

%Desired position and attitude
Fstar = 0; %desired thrust
euler_star = [0 0 -10]'.*pi/180;  %desired euler angle
euler_dstar = [0 0 0]';
tau_star = [0 0 0]';  %desired torque
w_star = [0 0 0]';
w_dstar = [0 0 0]';

%Body Frame
w = [0 0 0]';          % angular velocity in body frame
w_r = [0 0 0 0]';        % individual rotor speed
w_dot = [0 0 0]';
% Inertial frame
%pos = zeros(3,1000);        
pos(:,1)=[1 0 0];   % initial location
posdot = [0 0 0]';
posddot = [0 0 0]';
euler = [-10 5 10]'.*pi/180;
euler_dot  = [0 0 0]';

%Simulation

%R = [cos(euler(3,k))*cos(euler(2,k)) cos(euler(3,k))*sin(euler(2,k))*sin(euler(1,k))-sin(euler(3,k))*cos(euler(1,k)) cos(euler(3,k))*sin(euler(2,k))*cos(euler(1,k))+sin(euler(3,k))*sin(euler(1,k));
%     sin(euler(3,k))*cos(euler(2,k)) sin(euler(3,k))*sin(euler(2,k))*sin(euler(1,k))+cos(euler(3,k))*cos(euler(1,k)) sin(euler(3,k))*sin(euler(2,k))*cos(euler(1,k))-cos(euler(3,k))*sin(euler(1,k));
%     -sin(euler(2,k))                cos(euler(2,k))*sin(euler(1,k))                                                 cos(euler(2,k))*cos(euler(1,k))]
dt = 0.01;
for k = 1:2000
    
    ref(:,k) = [cos(k*dt);sin(k*dt);0.1*k*dt];  %fixed frame
    ref_dot(:,k) = [-sin(k*dt);cos(k*dt);0.1];  %fixed frame
    ref_ddot(:,k) = [-cos(k*dt);-sin(k*dt);0];  %fixed frame
    
    posddot(:,k) = ref_ddot(:,k) + Kv.*(ref_dot(:,k) - posdot(:,k)) + Kp.*(ref(:,k) - pos(:,k)); %commanded acc.
    Fstar(k) = (M*g + M*posddot(3,k))/(cos(euler(2,k))*cos(euler(1,k)));  % U1
    %desired eular angles
    euler_star(1,k+1) = asin((M*posddot(1,k)*sin(euler_star(3,k)) - M*posddot(2,k)*cos(euler_star(3,k)))/Fstar(k));
    euler_star(2,k+1) = asin((M*posddot(1,k)*cos(euler_star(3,k)) + M*posddot(2,k)*sin(euler_star(3,k)))/Fstar(k)*cos(euler_star(1,k)));
    euler_star(3,k+1) = euler_star(3,k);
    
    posdot(:,k+1) =  posdot(:,k) + dt*posddot(:,k);
    pos(:,k+1) = pos(:,k) + dt*posdot(:,k);
    
    euler_dot(:,k) = euler_dstar + Kpe.*(euler_star(:,k) - euler(:,k));
    euler(:,k+1) = euler(:,k) + dt*euler_dot(:,k);
    Retw = [1 0 -sin(euler(2,k));
            0 cos(euler(1,k)) cos(euler(2,k))*sin(euler(1,k)); 
            0 -sin(euler(1,k)) cos(euler(2,k))*cos(euler(1,k))];
    w(:,k) = Retw * euler_dot(:,k);
    
    w_dot(:,k) = w_dstar + Kpw.*(w_star - w(:,k));
                                 
    tau_star(:,k) = [Ixx*w_dot(1,k) + (Izz-Iyy)*w(2,k)*w(3,k); 
                     Iyy*w_dot(2,k) + (Ixx-Izz)*w(1,k)*w(3,k); 
                     Izz*w_dot(3,k) + (Iyy-Ixx)*w(1,k)*w(2,k)];
    w_r(:,k) = sqrt(inv([k k k k;
                         L*k 0 -L*k 0;
                         0 L*k 0 -L*k;
                         b -b b -b] ) * [Fstar(k);tau_star(:,k)]);
                     
    error_x(:,k) = ref(:,k) - pos(:,k);
end
t = 0:0.01:20;
figure(1)
plot(t,euler(3,:),t,euler(2,:),t,euler(1,:)),legend('\psi','\theta','\phi')

figure(2)
plot(t,pos(1,:),t,pos(2,:),t,pos(3,:)),legend('x','y','z')

t = linspace(0,20,2000);
%numel(t)
figure(3)
plot(t,tau_star(1,:),t,tau_star(2,:),t,tau_star(3,:));
legend("u2","u3","u4")
xlabel("time (s)")
ylabel("torque control input")
%ylim([-0.05 0.15])

t = linspace(0,20,2000);
%numel(t)
figure(4)
plot(t,Fstar);
legend("u1")
xlabel("time (s)")
ylabel("thrust control input")
%ylim([4.5 6])

figure(5)
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

figure(6)
plot(t,error_x(1,:),t,error_x(2,:),t,error_x(3,:)),legend('ex','ey','ez')
xlabel("time")
ylabel("error (m)")


% Z = stepinfo(pos(3,:),t,ref(3))
% Y = stepinfo(pos(2,:),t,ref(2))
% X = stepinfo(pos(1,:),t, ref(1))
