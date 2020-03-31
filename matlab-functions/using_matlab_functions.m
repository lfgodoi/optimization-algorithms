%%

% TITLE: 
   % MATLAB Optimization Functions

% DESCRIPTION:
   % Using MATLAB optimization functions.

% VERSION: 
   % Author: Leonardo Godói (eng.leonardogodoi@gmail.com)
   % Creation date: 22-April-2019

% REVISION HISTORY:
   % V1.0 | 22-April-2019 | Leonardo Godói | Creation

% ---------------------------------------------------------------------

% Settings
clear;
clc;

% Using 'fminsearch' and 'fminunc' for unconstrained minimization
options = optimoptions(@fminunc,'Display','off');
x = [0;0;0;0];
f1 = @(x) 10*(x(2)-x(1)^2)^2 + (1-x(1))^2;
x0 = [1;2];
fprintf("Minimization by 'fminsearch'\n");
x = fminsearch(f1,x0)
fprintf("\nMinimization by 'fminunc'\n");
x = fminunc(f1,x0,options)

% Using 'fmincon' for constrained minimization
options = optimoptions(@fmincon,'Display','off');
x0 = [1;2];
fprintf("\nMinimization by 'fmincon'\n");
x = fmincon(@(x)f2(x),x0,[],[],[],[],[],[],@(x)nonlcon(x),options)

% Using 'linProg' for linear optimization
A = [1  1  1  0  0  0  0  0  0;
     0  0  0  1  1  1  0  0  0;
     0  0  0  0  0  0  1  1  1;
     -1  0  0  -1  0  0  -1  0  0;
     0  -1  0  0  -1  0  0  -1  0;
     0  0  -1  0  0  -1  0  0  -1];
b = [2000;
     3000;
     1500;
     -2000;
     -2000;
     -1000];
lb = [0 0 0 0 0 0 0 0 0];
f = [25 20 30 30 25 25 20 15 23];
options = optimoptions('linProg','Display','off');
fprintf("\nLinear Optimization by 'linProg'\n");
[x,z] = linprog(f,A,b,[],[],lb,[],options)

% Objective function definition
function f = f2(x)
    f =  (x(1)-2)^4+(x(1)-2*x(2))^2;
end

% Constraints argumento definition
function [c,ceq] = nonlcon(x)
    ceq = 0;
    c = x(1)^2-x(2);
end

% ---------------------------------------------------------------------