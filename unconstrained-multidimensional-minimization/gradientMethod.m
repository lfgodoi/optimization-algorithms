%%

% TITLE: 
   % Gradient Method

% DESCRIPTION:
   % Function for unconstrained multidimensional minimization by 
   % Gradient Method.

% VERSION: 
   % Author: Leonardo Godói (eng.leonardogodoi@gmail.com)
   % Creation date: 27-March-2019

% REVISION HISTORY:
   % V1.0 | 27-March-2019 | Leonardo Godói | Creation

%% --------------------------------------------------------------

function x = gradientMethod(f,x0,ep,watch)
% 1st argument: (f) Objective function
% 2nd argument: (x0) Initial point
% 3rd argument: (ep) Tolerance for stopping criterion
% 4th argument: (watch) Flag for enabling report

    % Initializing variables
    cp = ep + 1;                % Stopping criterion value
    x = symvar(f);              % Objective function symbolic variable
    aux = sym('aux');           % Auxiliar symbolic variable for alpha
    grad(x) = gradient(f);      % Objective function gradient
    [~,m] = size(symvar(f));    % Output vector size          
    x = x0;                     % Initial point vector
    i = 0;                      % Iteration counter
    max = 1000;                 % Maximum number of iterations
    x_aux = num2cell(x);        % Auxiliar variable
    g = grad(x_aux{:});         % Initial value of gradient
    
    % Iterative sequence
    % Stopping criterion #1: Gradient norm
    % Stopping criterion #2: Number of iterations
    while cp > ep & i < max
        
        % Saving the previous point
        x_anterior = x;
        g_anterior = g;
        
        % Computing the new point
        x = x_anterior + aux*(-g_anterior);
        x_aux = num2cell(x);
        f_aux(aux) = f(x_aux{:});
        alpha = dichotomousSearch(f_aux,0,2,0.00001,false);
        
        % Unidimensional search step
        x = double(x_anterior + alpha*(-g_anterior));
        x_aux = num2cell(x);  
        g = double(grad(x_aux{:}));
        
        % Computing the value of stopping criterion
        cp = double(abs(norm(g)));
        
        % Printing current result
        if watch == true
            result = ['Iteration ',num2str(i),': ',num2str(x)];
            disp(result);
        end
        
        % Incrementing the counter
        i = i + 1;
        
    end

end

%% --------------------------------------------------------------