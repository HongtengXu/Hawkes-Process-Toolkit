%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Test different Hawkes process simulation methods.
% Given synthetic data generated via different simulation methods, we
% estimate the parameters via MLE, and plot the relative estimation errors
% w.r.t. the number of training sequences.
%
% Provider:
% Hongteng Xu @ Georgia Tech
% June 14, 2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
addpath('Simulation')
addpath('Learning')

options.N = 50; % the number of sequences
options.Nmax = 100; % the maximum number of events per sequence
options.Tmax = 100; % the maximum size of time window
options.tstep = 0.2;% the step length for computing sup intensity
options.M = 50; % the number of steps
options.GenerationNum = 5; % the number of generations
D = 4; % the dimension of Hawkes processes
nTest = 5;
nSeg = 5;
nNum = options.N/nSeg;

disp('Fast simulation of Hawkes processes with exponential kernel')
para1.mu = rand(D,1)/D;
para1.A = rand(D, D);
para1.A = 0.25 * para1.A./max(abs(eig(para1.A)));
para1.A = reshape(para1.A, [D, 1, D]);
para1.w = 1;
Seqs1 = SimulationFast_Thinning_ExpHP(para1, options);

disp('Thinning-based simulation of Hawkes processes with arbitrary kernel')
para2 = para1;
para2.kernel = 'exp';
para2.landmark = 0;
Seqs2 = Simulation_Thinning_HP(para2, options);

disp('Approximate simulation of Hawkes processes via branching process')
para3 = para2;
Seqs3 = Simulation_Branch_HP(para3, options);


disp('Learning Hawkes processes from synthetic data')
alg.LowRank = 0;
alg.Sparse = 0;
alg.GroupSparse = 0;
alg.outer = 4;
alg.rho = 0.1;
alg.inner = 4;
alg.thres = 1e-5;
alg.Tmax = [];



disp('Evaluation of quality of synthetic data.')

Err = zeros(3,nSeg,nTest);
for n = 1:nTest
    for i = 1:nSeg
        % initialize
        model.A = rand(D,1,D)./(D^2);
        model.mu = rand(D,1)./D;
        model.kernel = 'exp';
        model.w = 1;
        model.landmark = 0;

        model1 = model;
        model2 = model;
        model3 = model;
        model1 = Learning_MLE_Basis( Seqs1(1:i*nNum), model1, alg );
        model2 = Learning_MLE_Basis( Seqs2(1:i*nNum), model2, alg );
        model3 = Learning_MLE_Basis( Seqs3(1:i*nNum), model3, alg );

        Err(1,i,n) = norm([model1.mu; model1.A(:)] - [para1.mu; para1.A(:)])/...
            norm([para1.mu; para1.A(:)]);
        Err(2,i,n) = norm([model2.mu; model2.A(:)] - [para2.mu; para2.A(:)])/...
            norm([para2.mu; para2.A(:)]);
        Err(3,i,n) = norm([model3.mu; model3.A(:)] - [para3.mu; para3.A(:)])/...
            norm([para3.mu; para3.A(:)]);
    end
end

Error = mean(Err,3);
Std = std(Err, 0, 3);
figure
hold on
for i = 1:3
    errorbar(nNum:nNum:options.N, Error(i,:), Std(i,:), 'o-');
end
hold off
axis tight
xlabel('The number of training sequences');
ylabel('Relative estimation error')
legend('FastThinning', 'Thinning', 'Branching')
title('Learning results based on different simulation methods')
