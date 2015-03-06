function P = GeneratePoints()

alpha = pi*10.0/180.0 % rotation about x-axis
beta = 0.0            % rotation about y-axis
gamma = 0.0           % rotation about z-axis

offset = [0.0, 0.0, 0.0];

N = 200;
R = 50.0;
sigma = 2.0;


Rx = [1.0, 0.0, 0.0;
      0.0, cos(alpha), -sin(alpha);
      0.0, sin(alpha), cos(alpha)]; 

Ry = [cos(beta), 0.0, sin(beta);
      0.0, 1.0, 0.0;
      -sin(beta), 0.0, cos(beta)]; 
    
Rz = [cos(gamma), -sin(gamma), 0.0;
      sin(gamma), cos(gamma), 0.0;
      0.0, 0.0, 1.0];

theta = random('unif', 0, 2*pi, N, 1);
ERR = random('norm', 0, sigma, N, 3);

P = [R*cos(theta), R*sin(theta), zeros(N,1)];
P = (Rz*Ry*Rx*P')' + ones(N,1) * offset
P = P + ERR;

csvwrite('points.csv', P)
