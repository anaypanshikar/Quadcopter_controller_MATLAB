%%PHYSICAL QUANTITIES

global g m L Ixx Iyy Izz k b          % physical constants  
g = 9.81; 
m = 0.468; 
L = 0.225; 
Ixx = 4.856*10^-3;
Iyy = 4.856*10^-3;
Izz = 8.801*10^-3; 
k = 10^-6; 
b = 10^-7;
I = [Ixx, 0, 0; 0, Iyy, 0; 0, 0, Izz];

syms x xdot                        % Inertial frame
x = [0; 0; 0];
xdot = [0; 0; 0];

syms theta thetadot                % Inertial frame
theta = [0; 0; 0];
thetadot = [0; 0; 0];

 
%syms w1 w2 w3 w4
%inputs = [w1^2, w2^2, w3^2, w4^2];

%%SIMULATION

start_time = 0;
end_time = 20;
dt = 0.05;
times = start_time:dt:end_time;
N = numel(times);

figure
%hold on
for t = times
    % Take input from our controller.
    i = input(t);

    omega = thetadot_2_omega(theta, thetadot);

    % Compute linear and angular accelerations.
    a = acceleration(i, theta, m, g, k);
    omegadot = angular_acceleration(i, omega, I, L, b, k);
    omega = omega + dt * omegadot;
    thetadot = omega_2_thetadot(theta, omega); 
    theta = theta + dt * thetadot;
    xdot = xdot + dt * a;
    x = x + dt * xdot;
    %scatter(t, x(1), 15, 'filled', 'red') 
    scatter3(x(1), x(2), x(3), 15, 'filled', 'red')
    ylim([-20,20])
    xlim([-20,20])
    zlim([0,10])
    pause(0.00001)
    hold on
end

%%FUNCTIONS

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
    r = [1, sin(phi)*tan(th), cos(phi)*tan(th); 0, cos(phi), -sin(phi); 0, sin(phi)/cos(th), cos(phi)/cos(th)];
    omega = r\thetadot;
end

function thetadot = omega_2_thetadot(theta, omega)
    phi = theta(1,1);
    th = theta(2,1);
    psi = theta(3,1);

%transformation matrix for body to inertial - angular quantities
    r = [1, sin(phi)*tan(th), cos(phi)*tan(th); 0, cos(phi), -sin(phi); 0, sin(phi)/cos(th), cos(phi)/cos(th)];
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
global m g k
    if t<=7
        w1 = 1075;
        w2 = 1075;
        w3 = 1075;
        w4 = 1075;
    end
    if t>7
        w1 = m*g/4/k + 100*sin(pi*t);
        w2 = m*g/4/k + 100*sin(pi*t);
        w3 = m*g/4/k + 100*sin(pi*t);
        w4 = m*g/4/k + 100*sin(pi*t);
    end
    i = [w1^2, w2^2, w3^2, w4^2];
end
