%%

% TITLE: 
   % Penalty Method

% DESCRIPTION:
   % Function for constrained minimization by Penalty Method.

% VERSION: 
   % Author: Leonardo Godói (eng.leonardogodoi@gmail.com)
   % Creation date: 10-April-2019

% REVISION HISTORY:
   % V1.0 | 10-April-2019 | Leonardo Godói | Creation

%% --------------------------------------------------------------

function x = penaltyMethod(f,c_eq,c_ineq,x0,ep,watch)
% 1st argument: (f) Objective function
% 2nd argument: (c_eq) Equality constraints (c_eq = 0)
% 3rd argument: (c_ineq) Inequality constraints (c_ineq <= 0)
% 4th argument: (x0) initial point
% 5th argument: (ep) Tolerance for stopping criterion
% 6th argument: (watch) Flag for enabling report

    % Initializing the variables
    x = symvar(f);                         % Objective function variable
    [linhas,cols] = size(formula(c_ineq)); % Number of inequality constraints (rows)  
    v = zeros(linhas,1);                   % Activation vector (1 or 0) for the constraints  
    x = x0;                                % Initial point                              
    c = 40;                                % Initial penalty parameter
    alpha = 4;                             % Updating factor of penalty parameter
    cp = ep + 1;                           % Comparison value for stopping criterion
    i = 0;                                 % Iteration counter
    max = 100;                             % Iteration limiter
    
    % Iterative sequence
    % Stopping criterion #1: x variation norm
    % Stopping criterion #2: Number of iterations
    while cp > ep & i < max
        
        % Saving the previous point
        x_anterior = x;
        
        % Auxiliar variable
        x_aux = num2cell(x);
        
        % Inequality constraints according to the current values of x
        H = c_ineq(x_aux{:});
        
        % Verifying the inequality constraints
        % Activated (v = 1): Not satisfied.
        % Deactivated (v = 0): Satisfied.
        for k = 1:linhas
            if H(k) > 0
                v(k) = 1;
            else
                v(k) = 0;
            end
        end
        
        % Building the penalized function
        phi = f + c * sum(c_eq.*c_eq) + c * sum(v.*(c_ineq.*c_ineq));
        
        % Unconstrained minimization step
        x = newtonMethod(phi,x,ep);
        
        % Updating the penalty parameter
        c = c * alpha;
        
        % Incrementing the counter
        i = i + 1;
        
        % Updating the value of stopping criterion
        cp = norm(x - x_anterior);
        
        % Printing the iteration result
        if watch == true
            result = ['Iteration ',num2str(i),': '];
            disp(result);
            disp(x);
        end
        
    end

end

%% --------------------------------------------------------------