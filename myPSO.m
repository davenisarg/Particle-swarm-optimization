%%%%%%%%%% Particle Swarm optimization %%%%%%%%%%%%%%%%%
%%%%%%%%%% Author: Nisarg Dave %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Run Benchmark_PSO file %%%%%%%%%%%%%%%%%%%%%
%%% In Benchmark_PSO file uncomment the function that you wanna use %%%

function [xBest, fBest] = myPSO(fitnessfunc, xrange , InitialGuess_x,currentdirectionWeight,globalBestWeight,localBestWeight,total_population,maximum_iterations,Tolerance_for_Func,Tolerance_for_X,optimization_type,Startfrom_Guess,Weight_x0)

%%%%%% Must require function arguments are fitnessfunc and xrange %%%%%%%%
%%%% All other arguments are optional %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%getting the lower limit of xrange
LowerLimit_x = xrange(:,1);

%getting the upper limit of xrange
UpperLimit_x = xrange(:,2);

% calculate the dimension of datarange using the xrange limits
[n, m] = size(LowerLimit_x);

% Setting up InitialGuess_x 
if ~exist('InitialGuess_x','var')
    %InitialGuess_x = [];   % no direction
    %InitialGuess_x = [0.1*ones(n,1)];
    InitialGuess_x = 0.1*LowerLimit_x + 0.1*UpperLimit_x; %desired direction near global minima.
end

% setting up the current search direction weight
if ~exist('currentdirectionWeight','var')
    currentdirectionWeight = 0.9;
end

% setting up the Global best search direction weight
if ~exist('globalBestWeight','var')
    globalBestWeight = 0.9; 
end

% setting up the loacal best search direction weight
% setting up the Global best search direction weight
if ~exist('localBestWeight','var')
    localBestWeight = 0.6; 
end

% Initializing the total population as per data dimension
if ~exist('total_population','var')
    total_population = 1*n;
end

% Defining the maximum number of generations or iterations
if ~exist('maximum_iterations','var')
    maximum_iterations = 100;
end
% Tolerance factor for the function
% exit when variance in fitness function is < Tolerance factor of function
if ~exist('Tolerance_for_Func','var')
    Tolerance_for_Func = 1e-16; 
end

% Tolerance factor for location/distance of particle
% exit when variance in current state of particle is  < Tolerance for X values 
if ~exist('Tolerance_for_X','var')
    Tolerance_for_X = 1e-16; 
end

%setting up the type of the optimization we want to perform.
if ~exist('optimization_type','var')
    optimization_type = true;  %true for minimization, false for maximization
end

% Guessing the initial best point to start. Performs better random directed
% search if selected proper value.
if ~exist('Startfrom_Guess','var')
    Startfrom_Guess = true;
end

% Assigning weight of x0 for setting importance threshold from 0 to 1
if ~exist('Weight_x0','var')
    Weight_x0 = 0.95;
end

%finding the difference in limits (Maximum range of the search space)
Delta_X = UpperLimit_x - LowerLimit_x;

%printing result parameters
display = 'iter';
Print_iter = 1;

% Calculating constriction parameter K
phi = globalBestWeight + localBestWeight;
K = (2 / (abs(2 - phi - (sqrt(phi^2-4*phi)))));

% Deciding the optimization type
if optimization_type
    optimization_at = @min;
else
    optimization_at = @max;
end

% Initialize the population

% Creating random points as a starting coordinates for Particle swarms
m = total_population;  %population size
X1 = LowerLimit_x*ones(1,m) + ((UpperLimit_x-LowerLimit_x)*ones(1,m)).*rand(n,m);
X2 = LowerLimit_x*ones(1,m) + ((UpperLimit_x-LowerLimit_x)*ones(1,m)).*rand(n,m);


% Now The next step is to move the points to the InitialGuess using the
% provided Intitial weight parameter value
w = Weight_x0;
X0 = InitialGuess_x*ones(1,m);
X1 = w*X0 + (1-w)*X1;
X2 = w*X0 + (1-w)*X2;

% Setting positions
X = X1;
% setting Intitial "velocity" of the population
V = X2-X1;  

% If starting from random guess then setting up initial X and V
if Startfrom_Guess  % use for 
   X(:,1) = InitialGuess_x; 
   V(:,1) = zeros(size(InitialGuess_x));
end

%calculating updated lower and upper limit
X_Low = LowerLimit_x*ones(1,m);
X_Upp = UpperLimit_x*ones(1,m);
% calculating fitness for each particle in swarm
F = fitnessfunc(X);  

% best place of particle in search space
X_Best = X; 
% associated best fitness value for that particle at X_Best point
F_Best = F; 
% Global minima for given fitness function
[F_Global, I_Global] = optimization_at(F_Best);
% best point with respect to the global minimum.
X_Global = X(:, I_Global);

% Structured storage creation for the all variable parameters.
function Created_Structure = Creating_Structure(varargin)
    N_Inputs = length(varargin);
    for i=1:N_Inputs
        name = inputname(i);
        Created_Structure.(name) = varargin{i};
    end
end

% Creating data structure for saving the data for analysis
SavedData(maximum_iterations) = Creating_Structure(X, V, F, X_Best, F_Best, X_Global, F_Global, I_Global);

% Initializing variables into strcuture
SavedVariables.X_Global = zeros(n,maximum_iterations);
SavedVariables.F_Global = zeros(1,maximum_iterations);
SavedVariables.I_Global = zeros(1,maximum_iterations);
SavedVariables.X_Best_Var = zeros(n,maximum_iterations);
SavedVariables.F_Best_Var = zeros(1,maximum_iterations);
SavedVariables.X_Best_Mean = zeros(n,maximum_iterations);
SavedVariables.F_Best_Mean = zeros(1,maximum_iterations);
SavedVariables.X_Var = zeros(n,maximum_iterations);
SavedVariables.F_Var = zeros(1,maximum_iterations);
SavedVariables.X_Mean = zeros(n,maximum_iterations);
SavedVariables.F_Mean = zeros(1,maximum_iterations);
SavedVariables.iter = 1:maximum_iterations;

%%%%%%%%%%%%%%%%%  Particle Swarm Optimization %%%%%%%%%%%%%%%%

% Looping till maximum iterations
SavedVariables.exitFlag = 1;
for iter = 1:maximum_iterations
    
    % getting new particle points
    if iter > 1   % Updated particles
        r1 = rand(n,m);
        r2 = rand(n,m);
        V =  ...   % Update to the new V
             K*(currentdirectionWeight*V + ...    % Current search direction
             globalBestWeight*r1.*((X_Global*ones(1,m))-X) + ...  % Global direction
             localBestWeight*r2.*(X_Best-X));    % Local best direction
        %V =  ...   % Update to the new V
        %    currentdirectionWeight*V + ...    % Current search direction
        %    globalBestWeight*r1.*((X_Global*ones(1,m))-X) + ...  % Global direction
        %    localBestWeight*r2.*(X_Best-X);    % Local best direction
        X_New = X + V;  % Update position
        X = max(min(X_New, X_Upp), X_Low);   % Bounds to the search space
        F = fitnessfunc(X);   % Calculating new updated fitness
        F_Best_New = optimization_at(F_Best, F);   % Finding new best F
        idxUpdate = F_Best_New ~= F_Best;  % Reporting newly updated indices
        X_Best(:,idxUpdate) = X(:,idxUpdate);  % Reporting best postions
        F_Best = F_Best_New;
        [F_Global, I_Global] = optimization_at(F_Best); % Updating global and local best
        X_Global = X(:, I_Global); % best points ever
    end
    %%% Saving Structured Data
    SavedData(iter) = Creating_Structure(X, V, F, X_Best, F_Best, X_Global, F_Global, I_Global);
    SavedVariables.X_Global(:,iter) = X_Global;
    SavedVariables.F_Global(iter) = F_Global;
    SavedVariables.I_Global(iter) = I_Global;
    SavedVariables.X_Var(:,iter) = var(X, 0, 2);
    SavedVariables.X_Best_Var(:,iter) = var(X_Best, 0, 2);
    SavedVariables.X_Mean(:,iter) = mean(X, 2);
    SavedVariables.X_Best_Mean(:,iter) = mean(X_Best, 2);
    SavedVariables.F_Var(1,iter) = var(F);
    SavedVariables.F_Best_Var(1,iter) = var(F_Best);
    SavedVariables.F_Mean(1,iter) = mean(F);
    SavedVariables.F_Best_Mean(1,iter) = mean(F_Best);
    
    %%% Printing the final results
    xVar = norm(SavedVariables.X_Var(:,iter));
    if strcmp('iter',display)
        if mod(iter-1,Print_iter)==0
            fprintf('iter: %3d,  xBest: %7f,  fBest: %7f,   fVar: %7f,    xVar: %7f\n',...
                iter, SavedVariables.X_Global(iter),SavedVariables.F_Global(iter), SavedVariables.F_Var(1,iter),xVar);
        end
    end
    
    %%% Checking Convergence Criteria ( reaching the Tolerance)
    if SavedVariables.F_Var(1,iter) < Tolerance_for_Func
        SavedVariables.exitFlag = 0;
        SavedData = SavedData(1:iter);
        break
    elseif xVar < Tolerance_for_X
        SavedVariables.exitFlag = 2;
        SavedData = SavedData(1:iter);
        break
    end
end 

% Saving the final Xbest and fbest
xBest = SavedVariables.X_Global(:,end);
fBest = SavedVariables.F_Global(end);
SavedVariables.input = Creating_Structure(fitnessfunc, InitialGuess_x, LowerLimit_x, UpperLimit_x);  %Copy inputs
SavedVariables.fEvalCount = iter*m;
fprintf('Maximum Interations reached.\n');
end