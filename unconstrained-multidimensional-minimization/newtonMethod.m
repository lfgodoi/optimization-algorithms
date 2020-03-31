%%

% TITLE: 
   % Newton's Method

% DESCRIPTION:
   % Function for unconstrained multidimensional minimization by 
   % Newton's Method.

% VERSION: 
   % Author: Leonardo Godói (eng.leonardogodoi@gmail.com)
   % Creation date: 29-March-2019

% REVISION HISTORY:
   % V1.0 | 29-March-2019 | Leonardo Godói | Creation

%% --------------------------------------------------------------

function x = newtonMethod(f,x0,ep,watch)
% 1st argument: (f) Objective function
% 2nd argument: (x0) Initial point
% 3rd argument: (ep) Tolerance for stopping criterion
% 4th argument: (watch) Flag for enabling report

    % Initializing variables 
    cp = ep + 1;                   % Stopping criterion value
    x = symvar(f);                 % Objective function symbolic function
    h(x) = hessian(f);             % Objective function Hessian matrix  
    aux = sym('aux');              % Auxiliar symbolic variable for alpha         
    grad(x) = gradient(f);         % Objective function gradient       
    [~,m] = size(symvar(f));       % Output vector size          
    x = x0;                        % Initial point vector      
    i = 0;                         % Iteration counter     
    max = 100;                     % Maximum number of iterations 
    x_aux = num2cell(x);           % Auxiliar variables
    g = grad(x_aux{:});            % Initial value of gradient
    D = eye(m,m);                  % Initializing D as an identity matrix
    
    % Iterative sequence
    % Stopping criterion #1: Gradient norm
    % Stopping criterion #2: Number of iterations
    while i < max & cp > ep
        
        % Saving the previous point
        x_anterior = x;
        
        % Computing the descent direction
        d = -inv(h(x_aux{:}))*grad(x_aux{:});
        
        % Unidimensional search step
        x = x_anterior + aux*d;
        x_aux = num2cell(x);
        f_aux(aux) = f(x_aux{:});
        alpha = dichotomousSearch(f_aux,0,1,0.00001,false);
        
        % Computing the new point
        x = double(x_anterior + alpha*d);
        
        % Computing the value of stopping criterion
        x_aux = num2cell(x);
        cp = double(abs(norm(grad(x_aux{:}))));
        
        % Incrementing the counter
        i = i + 1;
        
    end

end

%% --------------------------------------------------------------