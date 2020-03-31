%%

% TITLE: 
   % Newton's Search

% DESCRIPTION:
   % Function for unconstrained unidimensional minimization by 
   % Newton's Search.

% VERSION: 
   % Author: Leonardo Godói (eng.leonardogodoi@gmail.com)
   % Creation date: 21-March-2019

% REVISION HISTORY:
   % V1.0 | 21-March-2019 | Leonardo Godói | Creation

%% --------------------------------------------------------------

function x = newtonSearch(f,ep,watch)
% 1st argument: (f) Objective function
% 2nd argument: (ep) Reference for stopping criterion
% 3rd argument: (watch) Flag for enabling report

    % Initializing variables 
    x = 0;        % Initial point of minimum
    cp = ep + 1;  % Stopping criterion value
    k = 1;        % Iteration counter

    % Computing first and second derivatives of objective function
    d1_f = diff(f,1);
    d2_f = diff(f,2);
    
    % Iterative sequence
    % Stopping criterion based on modulus of minimum derivative
    while cp >= ep
        
        % Updating x in i
        x_previous = x;
        
        % Updating x in i + 1
        x = double(x_previous - d1_f(x_previous)/d2_f(x_previous));
        
        % Analyzing the stopping criterion
        cp = abs(d1_f(x));
        
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