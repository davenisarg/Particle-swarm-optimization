function [y]=sumofpowers(x)
% this is an arbitrarily dimensioned function
% input: x - n x d input matrix (n d-dimensional vectors)

i = [1:size(x,2)]+1;
y = sum(bsxfun(@power,abs(x),i),2);