function [y]=ackley(x)
% this is an arbitrarily dimensioned function
% input: x - n x d input matrix (n d-dimensional vectors)
% global min is at x = [0,0,...,0]

a = 20; b = 0.2; c = 2*pi;
[n,d]=size(x);

y = -a*exp(-b*sqrt(1/d*sum(x.^2,2)))-exp(1/d*sum(cos(c*x),2))+a+exp(1);