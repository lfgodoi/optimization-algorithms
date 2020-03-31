%%

% TITLE: 
   % Dichotomous Search

% DESCRIPTION:
   % Function for unconstrained unidimensional minimization by 
   % Dichotomous Search.

% VERSION: 
   % Author: Leonardo Godói (eng.leonardogodoi@gmail.com)
   % Creation date: 18-March-2019

% REVISION HISTORY:
   % V1.0 | 18-March-2019 | Leonardo Godói | Creation

%% --------------------------------------------------------------

function x = dichotomousSearch(f,a,b,ep,watch)
% 1st argument: (f) Objective function
% 2nd argument: (a) Start of initial uncertainty range
% 3rd argument: (b) end of initial uncertainty range
% 4th argument: (ep) Distance between an intermediary and the midpoint
% 5th argument: (watch) Flag for enabling report

    % Initializing variables
    l = ep*2;        % Final length of uncertainty range
    k = 1;           % Iteration counter
    x = 0;           % Point of minimum
    cp = ep + 1;     % Stopping criterion value 
    max = 100;       % Maximum number of iterations

    % Analysis of uncertainty range length
    % Stopping criterion #1: Relative variation of objective function
    % Stopping criterion #2: Number of elapsed iterations
    while double(cp) >= double(ep) & k <= max
    
        % Updating the previous value
        x_previous = x;
        
        % Updating intermediary points
        lambda = ((a + b)/2) - ep;
        mi = ((a + b)/2) + ep;
    
        % Updating the uncertainty range
        if f(lambda) < f(mi)
            b = mi;
        else
            a = lambda;
        end
        
        % Computing the point of minimum
        x = (a + b)/2;
        
        % Analyzing the stopping criterion
        cp = abs(f(x) - f(x_previous))/(1 + abs(f(x_previous)));
        
        % Printing current result
        if watch == true
            result = ['Iteration ',num2str(k),': ',num2str(x)];
            disp(result);
        end
        
        % Incrementing counter
        k = k + 1;
        
    end
    
    % Plotting the objective function and indicating the point of minimum
    if watch == true
        clf;
        fplot(f,[-4,2],'Linewidth',2);
        hold on;
        fplot(x,f(x),'*','MarkerSize',10);
        text(x-0.5,f(x)+2,['x* = ',num2str(x)]);
    end
    
end

%% --------------------------------------------------------------

