global g m L Ixx Iyy Izz k b kzd kd kzp kp I        % physical constants  
g = 9.81; 
m = 0.468; 
L = 0.225; 
Ixx = 4.856*10^-3;
Iyy = 4.856*10^-3;
Izz = 8.801*10^-3; 
k = 10^-6; 
b = 10^-7;
I = [Ixx, 0, 0; 0, Iyy, 0; 0, 0, Izz];
kzd = 10; 
kd = 10; 
kzp = 10; 
kp = 10;

syms x xdot                        % Inertial frame
x = [0; 0; 2];
xdot = [0; 0; 0.5];

syms theta thetadot                % Inertial frame
theta = pi./180*[15; 15; 15];
deviation = 100;
thetadot = deg2rad(2 * deviation * rand(3, 1) - deviation);

 
%syms w1 w2 w3 w4
%inputs = [w1^2, w2^2, w3^2, w4^2];

start_time = 0;
end_time = 20;
dt = 0.1;
times = start_time:dt:end_time;
N = numel(times);

figure
hold on
for t = times
    % Take input from our controller.
    %i = input(t);
    [x, theta] = pd_controller(x, xdot, theta, thetadot);

    %omega = thetadot_2_omega(theta, thetadot);

    % Compute linear and angular accelerations.
    %a = acceleration(i, theta, m, g, k);
    %omegadot = angular_acceleration(i, omega, I, L, b, k);
    %omega = omega + dt * omegadot;
    %thetadot = omega_2_thetadot(theta, omega); 
    %theta = theta + dt * thetadot;
    %xdot = xdot + dt * a;
    %x = x + dt * xdot;
    %scatter(t, x(1), 15, 'filled', 'red') 
    %scatter3(x(1), x(2), x(3), 15, 'filled', 'red')
    scatter(t, 180./pi*theta(1), 15, 'filled', 'red')
    hold on
    scatter(t, 180./pi*theta(2), 15, 'filled', 'blue')
    scatter(t, 180./pi*theta(3), 15, 'filled', 'green')
    %scatter(t, x(3), 10, 'filled', 'red')
    ylim([-10,20])
    xlim([0,20])
    pause(0.00001)
    %hold on
end

function omegadot = angular_acceleration(inputs, omega, I, L, b, k)
    tau = torques(inputs, L, b, k);
    omegadot = inv(I) * (tau - cross(omega, I * omega));
end

function R = rotation(theta)
    phi = theta(1,1);
    th = theta(2,1);
    psi = theta(3,1);

%transformation matrix for body to inertial - linear quantities
    R = [cos(phi)*cos(psi) - cos(th)*sin(phi)*sin(psi), -cos(psi)*sin(phi) - sin(psi)*cos(phi)*cos(th), sin(th)*sin(psi); cos(th)*cos(psi)*sin(phi) + cos(phi)*sin(psi), cos(th)*cos(phi)*cos(psi) - sin(phi)*sin(psi), cos(psi)*sin(th); sin(phi)*sin(th), cos(phi)*sin(th), cos(th)];
end

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
    
function a = acceleration(inputs, theta, m, g, k)
    R = rotation(theta);
    T = R * thrust(inputs, k);
    a = [0; 0; -g] + 1 / m * T;
end

function T = thrust(inputs, k)
    % Inputs are values for ${\omega_i}^2$
    T = [0; 0; k * sum(inputs)];
end

function tau = torques(inputs, L, b, k)
    % Inputs are values for ${\omega_i}^2$
    tau = [
        L * k * (inputs(1) - inputs(3))
        L * k * (inputs(2) - inputs(4))
        b * (inputs(1) - inputs(2) + inputs(3) - inputs(4))
    ];
end

function i = input(t)
%global m g k
    if t<=7
        w1 = 1075;
        w2 = 1075;
        w3 = 1075;
        w4 = 1075;
    end
    if t>7
        w1 = 1075;
        w2 = 1075 - 100*sin(pi*t);
        w3 = 1075;
        w4 = 1075 + 100*sin(pi*t);
    end
    i = [w1^2, w2^2, w3^2, w4^2];
end

function [x, theta] = pd_controller(x, xdot, theta, thetadot)
    global m g k b Ixx Iyy Izz L kzd kp kd kzp I 
    dt = 0.1;
    
    % Controller gains, tuned by hand and intuition.
    %syms kzd = 2.5 kd = 1.75 kzp = 1.5 kp = 6
    
    phi = theta(1,1);
    th = theta(2,1);
    psi = theta(3,1);

    % Compute total thrust
    T = m*(g - kzd*xdot(3) + kzp*(1-x(3)))./(cos(th)*cos(phi));
    
    % Compute the torques
    tau1 = (-kd*thetadot(1) - kp*phi)*Ixx;
    tau2 = (-kd*thetadot(2) - kp*th)*Iyy;
    tau3 = (-kd*thetadot(3) - kp*psi)*Izz;
    
    tau = [tau1; tau2; tau3];
    
    % Solve for the inputs, $\gamma_i$
    gamma1 = T./(4*k) - tau2./(2*k*L) - tau3./(4*b);
    gamma2 = T./(4*k) - tau1./(2*k*L) + tau3./(4*b);
    gamma3 = T./(4*k) + tau2./(2*k*L) - tau3./(4*b);
    gamma4 = T./(4*k) + tau1./(2*k*L) + tau3./(4*b);
    
    i = [gamma1, gamma2, gamma3, gamma4];
    
    a = [0; 0; -g] + 1 ./ m * [0; 0; T];
    
    omega = thetadot_2_omega(theta, thetadot);
    omegadot = inv(I) * (tau - cross(omega, I * omega));
    omega = omega + dt * omegadot;
    thetadot = omega_2_thetadot(theta, omega);
    theta = theta + dt * thetadot;
    xdot = xdot + dt * a;
    x = x + dt * xdot;
    
end
