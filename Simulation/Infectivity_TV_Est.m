function At = Infectivity_TV_Est(t, model)

% Generate infectivity matrix A(t) in R^{U*U}
%
% T: the time interval
% t: current time
% A: U*L*U tensor, coefficients of A(t)
% option.sigma: bandwidth of Gaussian basis
% option.landmark: the location of Gaussian basis

U = size(model.A,3);
L = length(t);
At = zeros(U,L,U);

Basis = Kernel(t, model);

for u = 1:size(model.A,3)
    At(:,:,u) = model.A(:,:,u)*Basis';
end




