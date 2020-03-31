%%

% TITLE: 
   % Quadratic Regression Search

% DESCRIPTION:
   % Function for unconstrained unidimensional minimization by 
   % Quadratic Regression Search.

% VERSION: 
   % Author: Leonardo Godói (eng.leonardogodoi@gmail.com)
   % Creation date: 23-March-2019

% REVISION HISTORY:
   % V1.0 | 23-March-2019 | Leonardo Godói | Creation

%% --------------------------------------------------------------

function x = quadraticRegressionSearch(f,points,ep,watch)
% 1st argumento: (f) Objective function
% 2nd argumento: (points) Array containing three known points
% 3rd argumento: (ep) Reference for stopping criterion
% 4th argumento: (watch) Flag for enabling report

    % Initializing variables 
    x = 0;           % Point of minimum   
    cp = ep + 1;     % Stopping criterion value
    k = 1;           % Iteration counter
    
    % Reordering the points
    points = sort(points);
    points = [points(1) points(3) points(2)];
    
    % Computing the derivative of objective function
    d_f = derivada(f,1);
    
    % Iterative sequence
    % stopping criterion based on relative variation of function
    while cp >= ep    

        % Updating x in i
        x_previous = x;
        
        % Updating x in i + 1
        x = double((1/2)* ...
                (((points(2)^2-points(3)^2)*f(points(1))+ ...
                (points(3)^2-points(1)^2)*f(points(2))+ ...
                (points(1)^2-points(2)^2)*f(points(3))) ...
                / ...
                ((points(2)-points(3))*f(points(1))+ ...
                (points(3)-points(1))*f(points(2))+ ...
                (points(1)-points(2))*f(points(3)))));

        % Choosing new points for interpolation
        if x >= points(1) & x <= points(3) 
            if f(x) < f(points(3))
                points = [points(1) x points(3)];
            else
                points = [x points(3) points(2)];
            end
        else 
            if f(points(3)) < f(x)
                points = [points(1) points(3) x];
            else
                points = [points(3) x points(2)];
            end
        end
        
        % Reordering the points
        points = sort(points);
        points = [points(1) points(3) points(2)];
        
        % Analyzing stopping criterion
        cp = abs(f(x) - f(x_previous))/(1 + abs(f(x_previous)));
        
        % Apresentação do resultado parcial
        if watch == true
            result = ['Iteration ',num2str(k),': ',num2str(double(x))];
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
        text(x-3,f(x)+400,['x* = ',num2str(double(x))]);
    end
    
end

%% --------------------------------------------------------------