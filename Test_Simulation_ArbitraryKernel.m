%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Test different simulation methods of Hawkes processes with arbitrary
% kernels.
%
% Provider:
% Hongteng Xu @ Georgia Tech
% June 14, 2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
addpath('Simulation')
addpath('Learning')

options.N = 100; % the number of sequences
options.Nmax = 100; % the maximum number of events per sequence
options.Tmax = 30; % the maximum size of time window
options.tstep = 0.2;
options.M = 50;
options.GenerationNum = 5;
D = 4; % the dimension of Hawkes processes
nTest = 5;
nSeg = 5;
nNum = options.N/nSeg;



disp('Thinning-based simulation of Hawkes processes with arbitrary kernel')
para1.kernel = 'gauss';
para1.landmark = [0,3,6,9];
L = length(para1.landmark);
para1.mu = rand(D,1)/D;
para1.A = rand(D, D, L);
para1.A = 0.25 * para1.A./max(abs(eig(sum(para1.A,3))));
para1.A = reshape(para1.A, [D, L, D]);
para1.w = 2;
Seqs1 = Simulation_Thinning_HP(para1, options);

disp('Approximate simulation of Hawkes processes via branching process')
para2 = para1;
Seqs2 = Simulation_Branch_HP(para2, options);


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

Err = zeros(2,nSeg,nTest);
for n = 1:nTest
    for i = 1:nSeg
        % initialize
        model.A = rand(D,L,D)./(L*D^2);
        model.mu = rand(D,1)./D;
        model.kernel = 'gauss';
        model.w = 2;
        model.landmark = para1.landmark;

        model1 = model;
        model2 = model;
        model1 = Learning_MLE_Basis( Seqs1(1:i*nNum), model1, alg );
        model2 = Learning_MLE_Basis( Seqs2(1:i*nNum), model2, alg );

        Err(1,i,n) = norm([model1.mu; model1.A(:)] - [para1.mu; para1.A(:)])/...
            norm([para1.mu; para1.A(:)]);
        Err(2,i,n) = norm([model2.mu; model2.A(:)] - [para2.mu; para2.A(:)])/...
            norm([para2.mu; para2.A(:)]);
        
    end
end

Error = mean(Err,3);
Std = std(Err, 0, 3);
figure
hold on
for i = 1:2
    errorbar(nNum:nNum:options.N, Error(i,:), Std(i,:), 'o-');
end
hold off
axis tight
xlabel('The number of sequences');
ylabel('Estimation error')
legend('Thinning', 'Branching')
title('Learning results based on different simulation methods')
