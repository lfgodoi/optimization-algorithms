%%

% TITLE: 
   % Fibonacci Search

% DESCRIPTION:
   % Function for unconstrained unidimensional minimization by 
   % Fibonacci Search.

% VERSION: 
   % Author: Leonardo Godói (eng.leonardogodoi@gmail.com)
   % Creation date: 19-March-2019

% REVISION HISTORY:
   % V1.0 | 19-March-2019 | Leonardo Godói | Creation

%% --------------------------------------------------------------

function x = FibonacciSearch(f,a,b,n,watch)
% 1st argument: (f) Objective function
% 2nd argument: (a) Start of initial uncertainty range
% 3rd argument: (b) end of initial uncertainty range
% 4th argument: (n) Position in Fibonacci sequence
% 5th argument: (watch) Flag for enabling report

    % Iteration counter
    k = 1; 
    
    % Computing the value for position n in Fibonacci sequence
    Fn = zeros(n,1);
    Fn(1) = 1;
    Fn(2) = 1;
    for i = 3:n
            Fn(i) = Fn(i-1) + Fn(i-2); 
    end
    
    % Computing the final uncertainty value
    l = (b - a)/Fn(n);

    % Analysis of uncertainty range length
    % Stopping criterion #1: Relative variation of objective function
    % Stopping criterion #2: Avoiding negative indexes in sequence
    while b - a >= l & k < (n - 1)
    
        % Updating intermediary points
        lambda = a + (Fn(n-k-1)/Fn(n-k+1))*(b-a);
        mi = a + (Fn(n-k)/Fn(n-k+1))*(b-a);
    
        % Updating the uncertainty range
        if f(lambda) < f(mi)
            b = mi;
        else
            a = lambda;
        end
        
        % Computing the point of minimum
        x = (a + b)/2;
        
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


