function [y]=rosenbrock(x)
% this is an arbitrarily dimensioned function
% input: x - n x d input matrix (n d-dimensional vectors)
% global min is at x = [1,1,...,1]

y = sum(100*(x(:,2:end)-x(:,1:end-1).^2).^2+(x(:,1:end-1)-1).^2,2);