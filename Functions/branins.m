function [y]=branins(x)
% this is a 2-dimensional function
% input: x - n x 2 input matrix (n 2-dimensional vectors)
% 3 global mins
% (-pi, 12.275), (pi, 2.275), and (9.42478, 2.475)
% usually evaluated on x1 in [-5,10] and x2 in [0,15]

a = 1; b = 5.1/(4*pi^2); c = 5/pi; d = 6; e = 10; f = 1/(8*pi);
y = a*(x(:,2)-b*x(:,1).^2+c*x(:,1)-d).^2 + e*(1-f)*cos(x(:,1))+e;