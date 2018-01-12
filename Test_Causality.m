%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Learning Granger causality for Hawkes processes by MLE based method
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear


options.N = 2000; % the number of sequences
options.Nmax = 500; % the maximum number of events per sequence
options.Tmax = 200; % the maximum size of time window
options.tstep = 0.1;
options.dt = 0.1; % the length of each time step
options.M = 250; % the number of steps in the time interval for computing sup-intensity
options.GenerationNum = 100; % the number of generations for branch processing
D = 7; % the dimension of Hawkes processes



disp('Approximate simulation of Hawkes processes via branching process')
disp('Complicated gaussian kernel')
para1.kernel = 'gauss'; % the type of kernels per impact function
para1.w = 2; % the bandwidth of gaussian kernel
para1.landmark = 0:4:12; % the central locations of kernels
L = length(para1.landmark);

% initialize ground truth parameters
para1.mu = rand(D,1)/D;
para1.A = zeros(D, D, L);
for l = 1:L
    para1.A(:,:,l) = (0.5^l)*(0.5+ones(D));
end
mask = rand(D).*double(rand(D)>0.7);
para1.A = para1.A.*repmat(mask, [1,1,L]);
para1.A = 0.25*para1.A./max(abs(eig(sum(para1.A,3)))); % ensure the stationarity of Hawkes process
tmp = para1.A;
para1.A = reshape(para1.A, [D, L, D]);
for di = 1:D
    for dj = 1:D
        phi = tmp(di, dj, :);
        para1.A(dj, :, di) = phi(:);
    end
end
% two simulation methods
Seqs1 = Simulation_Branch_HP(para1, options);
% Seqs1 = Simulation_Thinning_HP(para1, options);


%%
% Visualize all impact functions and infectivity matrix
[A, Phi] = ImpactFunc( para1, options );



disp('Maximum likelihood estimation and basis representation') 
alg1.LowRank = 0; % without low-rank regularizer
alg1.Sparse = 1; % with sparse regularizer
alg1.alphaS = 1;
alg1.GroupSparse = 1; % with group-sparse regularizer
alg1.alphaGS = 100;
alg1.outer = 8;
alg1.rho = 0.1; % the initial parameter for ADMM
alg1.inner = 5;
alg1.thres = 1e-5;
alg1.Tmax = [];
alg1.storeErr = 0;
alg1.storeLL = 0;

model1 = Initialization_Basis(Seqs1); 
% learning the model by MLE
model1 = Learning_MLE_Basis( Seqs1, model1, alg1 ); 
[A1, Phi1] = ImpactFunc( model1, options );


% Visualize the infectivity matrix (the adjacent matrix of Granger causality graph)
figure
subplot(121)        
imagesc(A)
title('Ground truth of infectivity')
axis square
colorbar
subplot(122)        
imagesc(A1)
title('Estimated infectivity-MLE')
colorbar
axis square
