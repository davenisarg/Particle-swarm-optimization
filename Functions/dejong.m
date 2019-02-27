function [y,grad]=dejong(x)
% this is an arbitrarily dimensioned function
% input: x - n x d input matrix (n d-dimensional vectors)
% global min is at x = [0,0,...,0]

y = sum(x.^2,2);
grad = 2.*x;

