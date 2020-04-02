%%

% TITLE: 
   % Barrier Method

% DESCRIPTION:
   % Function for constrained minimization by Barrier Method.

% VERSION: 
   % Author: Leonardo Godói (eng.leonardogodoi@gmail.com)
   % Creation date: 12-April-2019

% REVISION HISTORY:
   % V1.0 | 12-April-2019 | Leonardo Godói | Creation

%% --------------------------------------------------------------

function x = barrierMethod(f,c_ineq,x0,ep,watch)
% 1st argument: (f) Objective function
% 2nd argument: (c_ineq) Inequality constraints (c_ineq <= 0)
% 3rd argument: (x0) initial point
% 4th argument: (ep) Tolerance for stopping criterion
% 5th argument: (watch) Flag for enabling report

    % Initializing variables
    cp = ep + 1;      % Reference for stopping criterion
    x = symvar(f);    % Objective function variable
    x = x0;           % Initial point
    c = 10;           % Barrier parameter
    alpha = 0.1;      % Updating factor of barrier parameter
    i = 0;            % Iteration counter
    max = 100;        % Iteration limiter
    c_inicial = c;    % Initial value of parameter
    ajuste = false;   % Flag to indicate if alpha has changed 
    
    % Iterative sequence
    % Stopping criterion #1: Variation of the point of minimum
    % Stopping criterion #2: number of iterations
    % NOTE: Criterion #1 is activated only if c did not restart
    while (cp > ep | ajuste == true) & i < max

        % Saving the previous point
        x_anterior = x;
        
        % Auxiliar variable for computing the function value at i
        x_aux = num2cell(x);
        
        % If there is any inequality constraint, compute P_des
        if isempty(c_ineq)
            P_ineq = 0;
        else
            [linhas,cols] = size(formula(c_ineq));   
            f_ineq = sym2cell(formula(c_ineq));
            for n = 1: linhas
                    P_ineq(n) = 1/f_ineq{n,1};
            end
        end
        
        % Defining the penalized problem    
        phi = f + c*(-sum(P_ineq));
        
        % Applying the unconstrained minimization
        x = newtonMethod(phi,x,ep);
        x_aux = num2cell(x);
        
        % Checking each constraint
        for k = 1:n
            rest = c_ineq(x_aux{:});
            % If any constraint has been violated... 
            % ... then adjust alpha and restart
            if rest(k) > 0 & ajuste == false
                x = x_anterior;
                c = c_inicial;
                alpha = alpha*2;
                ajuste = true;
                break;
            % If all constraints have been satisfied...
            % ... then update c normally
            else
                if k == n
                    c = c*alpha;
                    ajuste = false;
                end
            end
        end
        
        % If adjustment occcured, then restart the loop
        if ajuste ==  true
            continue;
        end
        
        % Updating the value of stopping criterion
        cp = norm(x - x_anterior);
        
        % Incrementing the counter
        i = i + 1;
        
        % Printing the iteration result
        if watch == true
            result = ['Iteration ',num2str(i),': '];
            disp(result);
            disp(x);
        end
        
    end
    
end

%% --------------------------------------------------------------
