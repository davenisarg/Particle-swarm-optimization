%%%%%%%%%% Artificial Bee Colony Optimization %%%%%%%%%%%%%%%%%
%%%%%%%%%% Author: Nisarg Dave %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Run Benchmark_ABC file %%%%%%%%%%%%%%%%%%%%%
%%% In Benchmark_ABC file uncomment the function that you wanna use %%%
function [bestbee, mincost] = myABC(fitnessfunc,range, nvar, max_iterations,acceleration,total_population)
    
%%%%%% Must require function arguments are fitnessfunc and range %%%%%%%%
%%%% All other arguments are optional %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % setting the maximum iterations or generations
    if ~exist('max_iterations','var')
        max_iterations = 100;
    end
    
    %total number of variables
    if ~exist('nvar','var')
        nvar = 2;
    end
    
    % setting up the total population
    if ~exist('total_population','var')
        total_population = 10;
    end
    
    %setting up initial acceleration
    if ~exist('acceleration','var')
        acceleration = 1;
    end
    
    varsize = [1 nvar]; 
    
    xmin = range(:,1);
    xmax = range(:,2);
    
    % Bounding limit
    lim = round(nvar*total_population);
    
    globest.cost = inf;
    
    %initializing location and cost parameters
    init.loc = [];
    init.cost = [];
    
    % creating bee swarm using initialization
    bee = repmat(init, total_population, 1);

    % Scout Bees initiate Food Source
    for i = 1:total_population
        bee(i).loc = unifrnd(xmin, xmax, varsize);
        bee(i).cost = fitnessfunc(bee(i).loc);

        if bee(i).cost < globest.cost
            globest = bee(i);
        end
    end
    
    % setting up iterations & respective best cost parameters
    C = zeros(total_population,1);
    bestcost = zeros(max_iterations,1);
    
    tempEBX = zeros(max_iterations,total_population);
    tempEBY = zeros(max_iterations,total_population);
    tempOBX = zeros(max_iterations,total_population);
    tempOBY = zeros(max_iterations,total_population);
    tempSBX = zeros(max_iterations,total_population);
    tempSBY = zeros(max_iterations,total_population);
    scout = [];
    
    % creating live visualization of Artificial Bee Colonies
    for iter = 1:max_iterations
        [X,Y] = meshgrid(xmin:0.1:xmax, xmin:0.1:xmax);
        
        % for Rastringin
        a = 10; nvar=2;
        contour(X,Y,a*nvar + X.^2 - a*cos(2*pi.*X) + Y.^2 - a*cos(2*pi.*Y), 10); hold on;
        
        % for Ackley
        %a=20; b=0.2; c= 2*pi;
        %contour(X,Y,-a*exp(-b*sqrt(1/nvar*(X.^2+Y.^2)))-exp(1/nvar*(cos(c*X)+cos(c*Y)))+a+exp(1),10); hold on;

        % for Dejong
        %contour(X,Y,X.^2 + Y.^2,10); hold on;
        
        % for sum of powers
        %contour(X,Y,((abs(X).^2) + (abs(Y).^3)),10); hold on;
        
        % for Gramacy & Lee
        %contour(X,Y,(((sin(10*pi*X) / (2*X)) + ((X-1).^4)) + ((sin(10*pi*Y) / (2*Y)) + ((Y-1).^4))) ,10); hold on;
        
        % for Rosenbrock
        %contour(X,Y,((100*(X(:,2:end)-X(:,1:end-1).^2).^2+(X(:,1:end-1)-1).^2) + (100*(Y(:,2:end)-Y(:,1:end-1).^2).^2+(Y(:,1:end-1)-1).^2)),10); hold on;
               
        scatter(0,0,35,'ok','filled');
        title('Optimization simulation of Benchmark Function');
        % Employed Bees
        for i = 1:total_population
            K = [1:i-1 i+1:total_population];
            k = K(randi([1 length(K)]));
            
            phi = acceleration*unifrnd(-1,1,varsize);
            
            newbee.loc = min(max(bee(i).loc + phi.*(bee(i).loc-bee(k).loc), xmin),xmax);
            newbee.cost = fitnessfunc(newbee.loc);
            
            if newbee.cost < bee(i).cost
                bee(i) = newbee;
            else
                C(i) = C(i) + 1;
            end
            tempEBX(iter,i) = bee(i).loc(1);
            tempEBY(iter,i) = bee(i).loc(2);
        end

        % Plot Employed Bees
        scatter(tempEBX(iter,:),tempEBY(iter,:),'xb');
        if (iter > 1)
            for i = 1:total_population              
                line([tempSBX(iter-1,i);tempEBX(iter,i)],[tempSBY(iter-1,i);tempEBY(iter,i)],'Color','b');
            end
        end
        frame(3*iter-2) = getframe(gcf); pause(0.0001);
        
        % Onlooker Bees
        F = zeros(total_population,1);
        for i = 1:total_population
            if (bee(i).cost >= 0)
                F(i) = 1/(1+bee(i).cost);
            else
                F(i) = 1+abs(bee(i).cost);                
            end
        end
        P = F/sum(F);

        for j = 1:total_population
            i=find(rand<=cumsum(P),1,'first');
            K = [1:i-1 i+1:total_population];
            k = K(randi([1 length(K)]));
            
            phi = acceleration*unifrnd(-1,1,varsize);
            
            newbee.loc = min(max(bee(j).loc + phi.*(bee(j).loc-bee(k).loc), xmin),xmax);
            newbee.cost = fitnessfunc(newbee.loc);
            
            if newbee.cost < bee(j).cost
                bee(j) = newbee;
            else
                C(j) = C(j) + 1;
            end
            tempOBX(iter,j) = bee(j).loc(1);
            tempOBY(iter,j) = bee(j).loc(2);
        end
        
        % Plot Onlooker Bees
        scatter(tempOBX(iter,:),tempOBY(iter,:),'xr');
        for i = 1:total_population
            line([tempEBX(iter,i);tempOBX(iter,i)],[tempEBY(iter,i);tempOBY(iter,i)],'Color','r', 'LineStyle', '--');
        end
        frame(3*iter-1) = getframe(gcf); pause(0.0001);        
        
        % Scout Bees
        for i = 1:total_population
            if C(i) >= lim
                scout = [scout ; i];
                bee(i).loc = unifrnd(xmin,xmax,nvar);
                bee(i).cost = fitnessfunc(bee(i).loc);
                C(i) = 0;
            end
            tempSBX(iter,i) = bee(i).loc(1);
            tempSBY(iter,i) = bee(i).loc(2);
        end
        
        % Plot Scout Bees
        scatter(tempSBX(iter,scout),tempSBY(iter,scout),'ok','filled');       
        for i = 1:length(scout)
            line([tempOBX(iter,scout(i));tempSBX(iter,scout(i))],[tempOBY(iter,scout(i));tempSBY(iter,scout(i))],'Color','k', 'LineWidth', 2);
        end
        scout = []; 
        frame(3*iter) = getframe(gcf); pause(0.0001);
        
        for i = 1:total_population
            if bee(i).cost < globest.cost
                globest = bee(i);
            end
        end
        
        bestcost(iter) = globest.cost;
        hold off;
        
        disp(['Iteration ' num2str(iter) ' | Minimum cost = ' num2str(bestcost(iter))] );
    end
    
    name = ['pop=' num2str(total_population)];
    
    clf
    for iter = 1:max_iterations
        semilogy(1:iter,bestcost(1:iter));
        axis([1 max_iterations 1e-7 1 ]);
        title(['Minimum Value Plot | Population = ' num2str(total_population) ' | MinVal = ' num2str(bestcost(iter))] );
        framec(iter) = getframe(gcf); pause(0.0001);        
    end
    
    % getting the best bee and mincost parameters
    name = ['Convergence  pop=' num2str(total_population)];
    bestbee = globest;
    mincost = bestcost;
end