% Let us continue with our regression problem from the linear and quadratic
% programming tutorial. 

%The problem boiled down to solving the problem minimize ||Ax-y|| for some
%suitable norm. Let us select the 2-norm, but this time solve the extension
%where we minimize the worst-case cost, under the assumption that the
%matrix A is uncertain,A=A+d, where ||d||2<1. This can be shown to be
%equivalent to minimizing ||Ax-y||2 + ||x||2.    

%This problem can easily be solved using YALMIP. Begin by defining the
%data.
clear
clc
x = [1 0 3 4 0 6]';
t = (0:0.02:2*pi)';
A = [sin(t) sin(2*t) sin(3*t) sin(4*t) sin(5*t) sin(6*t)];
e = (-4+8*rand(length(A),1));
y = A*x+e;

%As a first approach, we will do the modelling by hand, by adding second
%order cones using the low-level command cone. 

lambda = [0.01 0.05 0.1];
xhatArr=zeros(6,3);

xhat = sdpvar(6,1);
sdpvar u v
for i=1:3
    F = [cone(y-A*xhat,u), cone(lambda(i)*xhat,v)];
    solvesdp(F,u + v);
    xhatArr(:,i)=value(xhat);
end

error=zeros(3,1);
for i=1:3
    xhatArr(:,i)
    error(i,1)=norm(xhatArr(:,i)-x);
end
error
% F = [cone(y-A*xhat,u), cone(lambda*xhat,v)];
% solvesdp(F,u + v);

%By using the automatic modelling support in the nonlinear operator
%framework, we can alternatively write 

%solvesdp([],norm(y-A*xhat,2) + lambda*norm(xhat,2));

% YALMIP will automatically model this as a second order cone problem, and
% solve it as such if a second order cone programming solver is installed
% (SeDuMi, SDPT3 or Mosek). If no second order cone programming solver is
% found, YALMIP will convert the model to a semidefinite program and solve
% it using any installed semidefinite programming solver.    