%Data = load('fglobalupdatedrastrigin.mat');
Data1 = load('ABCle.mat');
figure
%X =(1:1:100)
%plot(X,Data.F_Global1,X,Data1.F_Global);
plot(Data1.F_Global1)
xlabel('Number of Iterations');
%ylabel("Fbest Convergence");
ylabel("F Global Convergence");
title("F Global best Vs Number of Iterations for Gramacy & Lee's Function");
%title("100 Fbest values for Gramacy & Lee's Benchmark function");