function y = sim_pnlm(model,u)
%% This function simulates the response of a parallel nolinearity model (polynom) to multi-inputs u

[~,nU] = size(u.dataSet);
nM = length(model);

if nU~=nM
    disp('The number of inputs must match the number of model pathways.')
    y = [];
    return;
else
    y = 0;
    for j = 1:nU
        %-- First, limit the input to the range of the polynomial model - to prevent large responses outside the range
        polyMin = model{j,1}.polyRange(1);
        polyMax = model{j,1}.polyRange(2);
        
        minExceedIndx = u.dataSet <= polyMin;
        u(minExceedIndx,j) = polyMin;

        maxExceedIndx = u.dataSet >= polyMax;
        u(maxExceedIndx,j) = polyMax;
        
        %-- Then, simulate the response
        y = y + nlsim(model{j,1},u(:,j));
    end
end

end