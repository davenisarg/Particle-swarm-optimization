clc; clear;

%InitialGuess_x = [0.1*ones(3,1)];  % No initial guess

% Minimize Ackley's function
a = 20; b = 0.2; c = 2*pi;
xrange = [-32.768*ones(3,1),32.768*ones(3,1)];
[n,d]=size(xrange);
fitnessfunc = @(x)(-a*exp(-b*sqrt(1/d*sum(x.^2,2)))-exp(1/d*sum(cos(c*x),2))+a+exp(1));

% Minimize Branin's function
%a = 1; b = 5.1/(4*pi^2); c = 5/pi; d = 6; e = 10; f = 1/(8*pi);
%xrange = [-5,10;0,15]
%[n,d]=size(xrange);
%fitnessfunc = @(x)(a*(x(:,2)-b*x(:,1).^2+c*x(:,1)-d).^2 + e*(1-f)*cos(x(:,1))+e);

% Minimize Dejong's function
%xrange = [-5.12*ones(10,1),5.12*ones(10,1)];
%fitnessfunc = @(x)(sum(x.^2,2));

% Minimize Rosenbrock's function
%xrange = [-5*ones(5,1),10*ones(5,1)] 
%fitnessfunc = @(x)(sum(100*(x(:,2:end)-x(:,1:end-1).^2).^2+(x(:,1:end-1)-1).^2,2));

% Minimize Rastringin's function
%xrange = [-5.12*ones(4,1),5.12*ones(4,1)]; 
%[n,d]=size(xrange);
%fitnessfunc = @(x)(10*d + sum(x.^2-10*cos(2*pi*x),2));

% Minimize SumofPowers' function
%xrange = [-ones(5,1),ones(5,1)]; %
%fitnessfunc = @(x)(sum(bsxfun(@power,abs(x),1),2));

% Minimize Gramacy & Lee function
%xrange = [0.5,2.5]; 
%fitnessfunc = @(x)((sin(10*pi*x) / (2*x)) + ((x-1).^4));   

%%%% Solving using myPSO

%Solving with default InitialGuess_X
[xBest, fBest] = myPSO(fitnessfunc,xrange);