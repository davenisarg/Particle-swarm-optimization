clear;
clf; 

% for Gramacy & Lee
%fitnessfunc = @(x)(((sin(10*pi*x(1)) / (2*x(1))) + ((x(1)-1).^4)) + ((sin(10*pi*x(2)) / (2*x(2))) + ((x(2)-1).^4))) ;   
%myABC(fitnessfunc,[0.5, 2.5])

% for Ackley
%a=20; b=0.2; c= 2*pi;nvar=2;
%fitnessfunc = @(x)(-a*exp(-b*sqrt(1/nvar*(x(1).^2 +x(2).^2)))-exp(1/nvar*(cos(c*x(1))+cos(c*x(2))))+a+exp(1))
%myABC(fitnessfunc,[-32.768, 32.768])

% for Rastringin
a = 10; nvar = 2;
fitnessfunc = @(x)(a*nvar + x(1).^2 - a*cos(2*pi.*x(1)) + x(2).^2 - a*cos(2*pi.*x(2)));   
myABC(fitnessfunc,[-5.12, 5.12])

% for Dejong
%fitnessfunc = @(x)(x(1).^2 + x(2).^2) ;
%myABC(fitnessfunc,[-65.536, 65.536])

% for Sum of powers
%fitnessfunc =  @(x)((abs(x(1)).^2) + (abs(x(2)).^3));
%myABC(fitnessfunc,[-1, 1])

% for Rosenbrock
%fitnessfunc = @(x)((100*(x(1)-x(1).^2).^2+(x(1)-1).^2)+(100*(x(2)-x(2).^2).^2+(x(2)-1).^2));
%myABC(fitnessfunc,[-5, 10])

