clear
clc
% sdpvar x1 x2 x3 x4
% L=set([4+x1+2*x2 x2+x3+2;
%        x2+x3+2 3*x2+x4]>0);
% see(L);
addpath(genpath('C:\Users\dohuan.ME197\Dropbox\Graduate Research(DB)\YALMIP'))
x1=sdpvar(1,1);
x2=sdpvar(1,1);

L=set([1-x1 x1+x2 x1;
       x1+x2 2-x2 0;
       x1 0 1+x2]>0)...
       +set((x1+x2)<1);
solvesdp(L,[]);