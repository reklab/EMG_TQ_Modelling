function [pnlm,yhat] = multi_input_pnlm(u,y,poly_orders,varargin)
%% This function identifies a parallel pathway of static nonlinearities,
%% where each nonlinearity is with respect to an input; i.e., this model:
%% y = f1(u1) + f2(u2) + ... + fn(un), where n is the number of inputs

%% Modifications: 
%-- EST: 2024-11-29: I realized that

options = {...
          };

if arg_parse(options,varargin)
    return
end

[~,n] = size(u);
Ts = y.domainIncr;

%% Initialize Ginputs as a cell array, where each element of the cell contains Chebychev basis expansions of input
Ginputs = cell(1,n); 

%% Before calculating Ginputs, normalize the inputs to [-1,1] interval
min_u = min(u.dataSet);
max_u = max(u.dataSet);
avg_u = (max_u + min_u)/2;
rng_u = max_u - min_u;
u_n = (u.dataSet - avg_u)*2./rng_u;

for j = 1:n
    Ginputs{1,j} = multi_tcheb(u_n(:,j),poly_orders(1,j));
end

%% We would like to solve the problem as a Least Sqaures in the form of A.X = B
%++ where A is the concatenation of basis functions of the inputs, excluding 
%++ the repeating order 0 (i.e., DC term) for inputs j = 2 to n
A = Ginputs{1,1};

for j = 2:n
    A = cat(2,A,Ginputs{1,j}(:,2:end));
end

%++ B is simply the output
B = y.dataSet;

%% Solve the least squares problem, deciding whether to demean the I/O or not as set of flag
[x,xSTD] = lscov(A,B);

%% 
yhat = A*x;
yhat = nldat(yhat,'domainIncr',Ts);

%% Produce output model as a cell array of nonlinearities
pnlm = cell(n,1);

%-- The first input is an expection
nlModel = polynom('polyType','tcheb');
n1 = 1;
n2 = poly_orders(1)+1;
set(nlModel,'polyCoef',x(n1:n2),'polyStd',xSTD(n1:n2),'polyRange',[min_u(1);max_u(1)]);
pnlm{1,1} = nlModel;

for j = 2:n
    nlModel = polynom('polyType','tcheb');
    n1 = n2+1;
    n2 = n1+poly_orders(j)-1;
    set(nlModel,'polyCoef',[0;x(n1:n2)],'polyStd',[0;xSTD(n1:n2)],'polyRange',[min_u(j);max_u(j)]);  %-- We add a 0 for the DC term.
    pnlm{j,1} = nlModel;
end

end