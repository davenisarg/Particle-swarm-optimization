function [y,grad]=rastrigin(x)
% this is an arbitrarily dimensioned function
% input: x - n x d input matrix (n d-dimensional vectors)
% global min is at x = [0,0,...,0]
[n,d]=size(x);

y = 10*d + sum(x.^2-10*cos(2*pi*x),2);
grad = 2.*x + 20*pi*sin(2*pi*x);