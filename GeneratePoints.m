function P = GeneratePoints()
%
%   P = GeneratePoints()
%
% This matlab script generates a random points on the 2D circle
% in the 3D space. The result is stored in points.csv
% To visualize the result, use the following Matlab command:
%
%   scatter3(P(:,1), P(:,2), P(:,3))
%


alpha = 180.0*pi/180.0 % rotation about x-axis
beta =  180.0*pi/180.0  % rotation about y-axis
gamma = 180.0*pi/180.0 % rotation about z-axis

offset = [100.0, 200.0, 300.0];

N = 15;
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

csvwrite('points_src.csv', P)

P = (Rz*Ry*Rx*P')' + ones(N,1) * offset
P = P + ERR;

csvwrite('points_dst.csv', P)
