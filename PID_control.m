% physical constants
global g m L Ixx Iyy Izz k b kzd kd kzp kp ki I kzi         
g = 9.81; 
m = 0.468; 
L = 0.225; 
Ixx = 4.856*10^-3;
Iyy = 4.856*10^-3;
Izz = 8.801*10^-3; 
k = 10^-6; 
b = 10^-7;
I = [Ixx, 0, 0; 0, Iyy, 0; 0, 0, Izz];

%Controller gains
kzd = 5; 
kd = 5; 
kzp = 10; 
kp = 7;
ki = 12;
kzi = 10;

syms x xdot                        % Inertial frame
x = [0; 0; 2];                     %some initial location which needs to be stabilized to (0,0,1)
xdot = [0; 0; 0.5];

syms theta thetadot                % Inertial frame
theta = pi./180*[-10; 10; 5];      % some initial euler angles, which need to be driven to zero
deviation = 100;
thetadot = deg2rad(2 * deviation * rand(3, 1) - deviation);

%Simulation
start_time = 0;
end_time = 20;
dt = 0.1;
times = start_time:dt:end_time;
N = numel(times);

figure
hold on
for t = times
    % Take input from our controller.
    [x, theta] = pd_controller(x, xdot, theta, thetadot, t);
    
    %plot for euler angles:
    %scatter(t, 180./pi*theta(1), 15, 'filled', 'red')
    %hold on
    %scatter(t, 180./pi*theta(2), 15, 'filled', 'blue')
    %scatter(t, 180./pi*theta(3), 15, 'filled', 'green')
    
    %plot for position:
    scatter(t, x(3), 15, 'filled', 'red')
    scatter(t, x(2), 15, 'filled', 'blue')
    scatter(t, x(1), 15, 'filled', 'green')
    legend('z','y','x')
    ylim([-0.5,2])                  %take [-2,3] for z coord plot
    xlim([0,10])                       
    pause(0.00001)
    %hold on
end

%Functions
function omega = thetadot_2_omega(theta, thetadot)
    phi = theta(1,1);
    th = theta(2,1);
    psi = theta(3,1);

%transformation matrix for body to inertial - angular quantities
    r = [1, sin(phi)*tan(th), cos(phi)*tan(th); 0, cos(phi), -sin(phi); 0, sin(phi)./cos(th), cos(phi)./cos(th)];
    omega = r\thetadot;
end

function thetadot = omega_2_thetadot(theta, omega)
    phi = theta(1,1);
    th = theta(2,1);
    psi = theta(3,1);

%transformation matrix for body to inertial - angular quantities
    r = [1, sin(phi)*tan(th), cos(phi)*tan(th); 0, cos(phi), -sin(phi); 0, sin(phi)./cos(th), cos(phi)./cos(th)];
    thetadot = r*omega;
end

function [x, theta] = pd_controller(x, xdot, theta, thetadot, t)
    global m g k b Ixx Iyy Izz L kzd kp kd ki kzp I kzi
    dt = 0.1;
    
    % Controller gains, tuned by hand and intuition.
    %syms kzd = 2.5 kd = 1.75 kzp = 1.5 kp = 6
    z = x(3);
    phi = theta(1,1);
    th = theta(2,1);
    psi = theta(3,1);

    %PID CONTROLLER FOR Z COORDINATE: STABILIZES TO LOCATION (0,0,1)
    ez = 0;
    t1 = 0;
    while t1 <= t
        ez = ez + (1-z)*dt;
        t1 = t1 + dt;
    end
    etz = - kzd*xdot(3) + kzp*(1-x(3)) + kzi*ez;
    T = m*(g + etz)./(cos(th)*cos(phi));
    
    %PID CONTROLLER FOR EULER ANGLES: DRIVES ALL TO ZERO
    e1 = 0;
    e2 = 0;
    e3 = 0;
    t1 = 0;
    while t1 <= t
        e1 = e1 + -phi*dt;
        e2 = e2 + -th*dt;
        e3 = e3 + -psi*dt;
        t1 = t1 + dt;
    end
    et1 = -kd*thetadot(1) - kp*phi + ki*e1;
    et2 = -kd*thetadot(2) - kp*th + ki*e2;
    et3 = -kd*thetadot(3) - kp*psi + ki*e3;
    tau1 = (et1)*Ixx;
    tau2 = (et2)*Iyy;
    tau3 = (et3)*Izz;
    
    tau = [tau1; tau2; tau3];
    
    % Solve for the inputs, $\gamma_i$
    gamma1 = T./(4*k) - tau2./(2*k*L) - tau3./(4*b);
    gamma2 = T./(4*k) - tau1./(2*k*L) + tau3./(4*b);
    gamma3 = T./(4*k) + tau2./(2*k*L) - tau3./(4*b);
    gamma4 = T./(4*k) + tau1./(2*k*L) + tau3./(4*b);
    
    %i = [gamma1, gamma2, gamma3, gamma4];
    
    a = [0; 0; -g] + 1 ./ m * [0; 0; T];
    
    omega = thetadot_2_omega(theta, thetadot);
    omegadot = inv(I) * (tau - cross(omega, I * omega));
    omega = omega + dt * omegadot;
    thetadot = omega_2_thetadot(theta, omega);
    theta = theta + dt * thetadot;
    xdot = xdot + dt * a;
    x = x + dt * xdot;
    
end
