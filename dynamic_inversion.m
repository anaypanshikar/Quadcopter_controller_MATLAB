close all
clear all 
clc 
% physical constants     
g = 9.81; 
M = 0.5; 
L = 0.25;
Ixx = 5*10^-3;
Iyy = 5*10^-3;
Izz = 10*10^-3; 
k = 3*10^-6; 
b = 10^-7;
I = diag([Ixx Iyy Izz]);

%constants
Kv = [0.63 ; 0.63; 10.5];  %for positions
Kp = [0.1225 ; 0.1225; 6.125];
Kpe = [6; 5.2; 7];      %for euler angles
Kpw = [4; 3; 4.9];       % for w

%Desired position and attitude
ref = [5 6 4]';
ref_dot = [0 0 0]';
ref_ddot = [0 0 0]';
Fstar = 0.5*9.81; %desired thrust
euler_star = [0 0 10]'.*pi/180;  %desired euler angle
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
pos(:,1)=[0 0 1];   % initial location
posdot = [0 0 0]';
posddot = [0 0 0]';
euler = [10 -10 0]'.*pi/180;
euler_dot  = [0 0 0]';

%Simulation

%R = [cos(euler(3,k))*cos(euler(2,k)) cos(euler(3,k))*sin(euler(2,k))*sin(euler(1,k))-sin(euler(3,k))*cos(euler(1,k)) cos(euler(3,k))*sin(euler(2,k))*cos(euler(1,k))+sin(euler(3,k))*sin(euler(1,k));
%     sin(euler(3,k))*cos(euler(2,k)) sin(euler(3,k))*sin(euler(2,k))*sin(euler(1,k))+cos(euler(3,k))*cos(euler(1,k)) sin(euler(3,k))*sin(euler(2,k))*cos(euler(1,k))-cos(euler(3,k))*sin(euler(1,k));
%     -sin(euler(2,k))                cos(euler(2,k))*sin(euler(1,k))                                                 cos(euler(2,k))*cos(euler(1,k))]
dt = 0.01;
for k = 1:2000
   
    posddot(:,k) = ref_ddot + Kv.*(ref_dot - posdot(:,k)) + Kp.*(ref - pos(:,k)); %commanded acc.
    Fstar(k+1) = (M*g + M*posddot(3,k))/(cos(euler(2,k))*cos(euler(1,k)));  % U1
    %desired euler angles
    euler_star(1,k+1) = asin((M*posddot(1,k)*sin(euler_star(3,k)) - M*posddot(2,k)*cos(euler_star(3,k)))/Fstar(k));
    euler_star(2,k+1) = asin((M*posddot(1,k)*cos(euler_star(3,k)) + M*posddot(2,k)*sin(euler_star(3,k)))/Fstar(k)*cos(euler_star(1,k)));
    euler_star(3,k+1) = euler_star(3,k);
    
    posdot(:,k+1) =  posdot(:,k) + dt*posddot(:,k);
    pos(:,k+1) = pos(:,k) + dt*posdot(:,k);
    
    euler_dot(:,k+1) = euler_dstar + Kpe.*(euler_star(:,k) - euler(:,k));         %%%%%
    euler(:,k+1) = euler(:,k) + dt*euler_dot(:,k);
    Retw = [1 0 -sin(euler(2,k));
            0 cos(euler(1,k)) cos(euler(2,k))*sin(euler(1,k)); 
            0 -sin(euler(1,k)) cos(euler(2,k))*cos(euler(1,k))];
    w(:,k) = Retw * euler_dot(:,k);
    
    w_dot(:,k+1) = w_dstar + Kpw.*(w_star - w(:,k));
                                 
    tau_star(:,k) = [Ixx*w_dot(1,k) + (Izz-Iyy)*w(2,k)*w(3,k); 
                     Iyy*w_dot(2,k) + (Ixx-Izz)*w(1,k)*w(3,k); 
                     Izz*w_dot(3,k) + (Iyy-Ixx)*w(1,k)*w(2,k)];
    w_r(:,k) = sqrt(inv([k k k k;
                         L*k 0 -L*k 0;
                         0 L*k 0 -L*k;
                         b -b b -b] ) * [Fstar(k);tau_star(:,k)]);
    error(:,k) = ref - pos(:,k); 
end
t = 0:0.01:20;
figure
plot(t,euler(1,1:2001)*180/pi,t,euler(2,1:2001)*180/pi,t,euler(3,1:2001)*180/pi),legend('\phi','\theta','\psi')
xlabel('time')
ylabel('euler angles (degrees)')
ylim([-15 15])

t = 0:0.01:10;
figure
plot(t,euler_dot(1,1:1001)*180/pi,t,euler_dot(2,1:1001)*180/pi,t,euler_dot(3,1:1001)*180/pi),legend('$\dot{\phi}$','$\dot{\theta}$','$\dot{\psi}$','interpreter','latex')
xlabel('time')
ylabel('angular velocities (deg/s)')

t = 0:0.01:20;
figure
plot(t,pos(1,:),t,pos(2,:),t,pos(3,:)),legend('x','y','z')
xlabel('time')
ylabel('position (m)')
ylim([0 8])

t = 0:0.01:20;
figure
plot(t,Fstar),legend('U1')
xlabel('time')
ylabel('Thrust (N)')

t = linspace(0,20,2000);
figure
plot(t,tau_star(1,:),t,tau_star(2,:),t,tau_star(3,:)),legend('U2',"U3","U4")
xlabel('time')
ylabel('Torques (Nm)')

t = linspace(0,20,2000);
figure
plot(t,error(1,:),t,error(2,:),t,error(3,:)),legend('ex',"ey","ez")
xlabel('time')
ylabel('Error (m)')

%Z = stepinfo(pos(3,:),t,ref(3))
%Y = stepinfo(pos(2,:),t,ref(2))
%X = stepinfo(pos(1,:),t, ref(1))
