%%

% TITLE: 
   % Quasi-Newton Method

% DESCRIPTION:
   % Function for unconstrained multidimensional minimization by 
   % Quasi-Newton Method.

% VERSION: 
   % Author: Leonardo Godói (eng.leonardogodoi@gmail.com)
   % Creation date: 30-March-2019

% REVISION HISTORY:
   % V1.0 | 30-March-2019 | Leonardo Godói | Creation

%% --------------------------------------------------------------

function x = quasiNewtonMethod(f,x0,ep)
% 1st argument: (f) Objective function
% 2nd argument: (x0) Initial point
% 3rd argument: (ep) Tolerance for stopping criterion

    % Initializing variables                                   
    cp = ep + 1;                    % Stopping criterion value                
    aux = sym('aux');               % Auxiliar symbolic variable for alpha         
    x = symvar(f);                  % Objective function symbolic variable
    grad(x) = gradient(f);          % Objective function gradient       
    [~,m] = size(symvar(f));        % Output vector size          
    x = x0;                         % Initial point vector      
    i = 0;                          % Iteration counter     
    max = 100;                      % maximum number of iterations 
    x_aux = num2cell(x);            % Auxiliar variables
    g = double(grad(x_aux{:}));     % Initial value of gradient
    D = eye(m,m);                   % Initializing D as identity

    % Iterative sequence
    % Stopping criterion #1: Gradient norm
    % Stopping criterion #2: Number of iterations
    while cp > ep & i < max

        % Saving the previous point
        x_anterior = x;
        g_anterior = g;

        % Computing the descent direction
        d = -D*g;

        % Unidimensional search step
        x = x_anterior + aux*d;
        x_aux = num2cell(x);
        f_aux(aux) = f(x_aux{:});
        alpha = dichotomousSearch(f_aux,0,1,0.00001,false);

        % Computing the new point
        x = double(x_anterior + alpha*d);
        x_aux = num2cell(x);
        g = double(grad(x_aux{:}));

        % Evaluating stopping criterion
        cp = abs(norm(g));

        % If not satisfied, then update D
        if cp > ep
            s = double(x - x_anterior);
            y = double(g - g_anterior);
            if s'*y <= 0
                D = eye(m,m);
            else
                D = D + ((s*s')/(s'*y)) - ((D*y*y'*D)/(y'*D*y));
            end
        end

        % Incrementing the counter
        i = i + 1;

    end

end

%% --------------------------------------------------------------
