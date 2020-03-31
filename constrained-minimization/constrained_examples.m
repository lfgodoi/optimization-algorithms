%%

% TITLE: 
   % Constrained Optimization Examples

% DESCRIPTION:
   % Applying all the constrained minimization methods implemented.

% VERSION: 
   % Author: Leonardo Godói (eng.leonardogodoi@gmail.com)
   % Creation date: 12-April-2019

% REVISION HISTORY:
   % V1.0 | 12-April-2019 | Leonardo Godói | Creation

%% --------------------------------------------------------------

% Penalty Method
clear;
clc;
ep = 0.00001;   
x = sym('x',[2 1]); 
f(x) = (1/2)*(x(1)^2+(x(2)^2/3));
h = symfun([0],x);
c = symfun([(x(1)+x(2)-1)],x); 
x0 = [0;0];
fprintf("Result by Penalty Method:\n")
penaltyMethod(f,c,h,x0,ep,false)

% Barrier Method
ep = 0.00001;   
x = sym('x',[2 1]); 
f(x) = (x(1)-2)^4+(x(1)-2*x(2))^2;
r_des = symfun([x(1)^2-x(2)], x);
x0 = [1;2];
fprintf("\nResult by Barrier Method:\n")
barrierMethod(f,r_des,x0,ep,false)

% Augmented Lagrangian Method
ep = 0.00001;   
x = sym('x',[2 1]); 
f(x) = (1/2)*(x(1)^2+(x(2)^2/3));
h = symfun([0],x);
c = symfun([(x(1)+x(2)-1)],x); 
x0 = [0;0];
fprintf("\nResult by Augmented Lagrangian Method:\n")
augmentedLagrangianMethod(f,c,h,x0,ep,false)

%% --------------------------------------------------------------