syms s kp1 
kp1 = 46.34;

%TF
G = tf([kp1],[1,0,0]);

margin(G)

step(G/(1+G))
xlim([0 50])

K = 1;
P = 0.923;
kp = 0.6*K;
ti = 0.5*P;
td = 0.125*P;

s = tf('s');
Gc = K*(1 + 1/(ti*s) + td*s);
Go = kp1 / s^2;
G1 = Gc*Go;
%margin(G1)
step(G1/(1+G1))

