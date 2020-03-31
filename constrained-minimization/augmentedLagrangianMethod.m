%%

% TITLE: 
   % Augmented Lagrangian Method

% DESCRIPTION:
   % Function for constrained minimization by Augmented Lagrangian Method.

% VERSION: 
   % Author: Leonardo Godói (eng.leonardogodoi@gmail.com)
   % Creation date: 08-April-2019

% REVISION HISTORY:
   % V1.0 | 08-April-2019 | Leonardo Godói | Creation

%% --------------------------------------------------------------

function x = augmentedLagrangianMethod(f,c,h,x0,ep,watch)
% 1st argument: (f) Objective function
% 2nd argument: (c) Equality constraints
% 3rd argument: (h) Equality constraints
% 4th argument: (x0) Initial point
% 5th argument: (ep) Tolerance for stopping criterion
% 6th argument: (watch) Flag for enabling report

    % Initializing variables
    cp = ep + 1;               % Reference for stopping criterion
    x = symvar(f);             % Objective function variable
    x = x0;                    % Initial point
    nu = zeros(size(c));       % Lagrange Multiplier for equality
    lambda = zeros(size(h));   % Lagrange Multiplier for inequality
    ri = ones(size(c));        % Penalty parameters for equality
    rj = ones(size(c));        % Penalty parameters for equality
    k = 0;                     % Iteration counter
    max = 100;                 % Maximum number of iterations
    tol = ep/100;              % Tolerance for updating Lagrange multipliers
    alpha = 2;                 % Penalty updating parameter
    
    % Iterative sequence
    % Stopping criterion #1: Variation of the point of minimum
    % Stopping criterion #2: number of iterations
    while cp > ep & k < max

        % Saving the previous point
        x_anterior = x;
        
        % Defining the penalized problem
        phi = f + sum(nu.*c) + sum(lambda.*h)...
              + (1/2)*sum(ri.*(c.*c)) + (1/2)*sum(rj.*(h.*h));
        
        % Applying the unconstrained minimization
        x = double(quasiNewtonMethod(phi,x,ep));
        
        % Auxiliar variable for computing the function value at i
        x_aux = num2cell(x);
        
        % Updating the multipliers
        c_aux = double(c(x_aux{:}));
        h_aux = double(h(x_aux{:}));
        for i = 1:size(c)
            % Equality constraints
            if abs(c_aux(i)) > tol
                nu(i) = nu(i) + ri(i)*c_aux(i);
                ri(i) = alpha*ri(i);
            end
        end
        for j = 1:size(h)
            % Inequality constraints
            if h_aux(j) > tol
                lambda(j) = lambda(j) + rj(j)*h_aux(j);
                rj(j) = alpha*rj(j);
            end
        end
        
        % Updating the value of stopping criterion
        cp = norm(x - x_anterior);
        
        % Incrementing the counter
        k = k + 1;
        
        % Printing the iteration result
        if watch == true
            result = ['Iteration ',num2str(k),': '];
            disp(result);
            disp(x);
        end
    
    end
    
end

%% --------------------------------------------------------------
