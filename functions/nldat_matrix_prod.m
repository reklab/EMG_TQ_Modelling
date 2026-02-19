function y = nldat_matrix_prod(u,M)
%% This function performs a matrix product of a multi-channel nldat signal with a matrix from right side.
%% Therefore, it creates a linear combination of nldat channels.
%% Inputs: 
%       + u is an nldat object of size [N,nInputs]
%       + M is a matrix of size [nInputs,nOutputs]

%% Output:
%       + y is an nldat object of size [N,nOutputs]

[~,nInputs] = size(u.dataSet);
[rM,~] = size(M);
Ts = u.domainIncr;

if rM < nInputs
    disp('The matrix M has less number of rows than the signals in u!')
    y = [];
    return;
elseif rM > nInputs
    disp('The matrix M has more number of rows than the signals in u!')
    y = [];
    return;
else
    y = nldat(u.dataSet * M,'domainIncr',Ts);
end
