%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compare batch opt of Hawkes with online/stochastic opt of Hawkes
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
% data simulation
options.N = 50; % the number of sequences
options.Nmax = 1000; % the maximum number of events per sequence
options.Tmax = 100; % the maximum size of time window
options.tstep = 0.2;% the step length for computing sup intensity
options.M = 50; % the number of steps
options.GenerationNum = 5; % the number of generations
D = 10; % the dimension of Hawkes processes

nTest = 5;
Nout = 40;

disp('Fast simulation of Hawkes processes with exponential kernel')
para.mu = rand(D,1)/D;
para.A = rand(D, D);
para.A = 0.65 * para.A./max(abs(eig(para.A)));
para.A = reshape(para.A, [D, 1, D]);
para.w = 1;
Seqs = SimulationFast_Thinning_ExpHP(para, options);

err1 = zeros(nTest, Nout);
err2 = zeros(nTest, Nout);

for n = 1:nTest

% initialize
model.A = rand(D,1,D)./(D^2);
model.mu = rand(D,1)./D;
model.kernel = 'exp';
model.w = 1;
model.landmark = 0;

disp('Learning HP by batch opt')
alg1.LowRank = 0;
alg1.Sparse = 0;
alg1.GroupSparse = 0;
alg1.outer = Nout;
alg1.rho = 0.1;
alg1.inner = 1;
alg1.thres = 1e-5;
alg1.Tmax = [];
alg1.storeErr = 1;
alg1.storeLL = 0;
alg1.truth = para;
model1 = Learning_MLE_Basis( Seqs, model, alg1 );

disp('Learning HP by stochastic opt')
alg2.LowRank = 0;
alg2.Sparse = 0;
alg2.GroupSparse = 0;
alg2.epoch = Nout;
alg2.rho = 0.1;
alg2.eventbatch = 20;
alg2.seqbatch = 10;
alg2.historyL = 20;
alg2.thres = 1e-5;
alg2.Tmax = [];
alg2.storeErr = 1;
alg2.storeLL = 0;
alg2.truth = para;
model2 = Learning_MLE_Basis_Stoc( Seqs, model, alg2 );

err1(n,:) = model1.err(:,3)';
err2(n,:) = model2.err(:,3)';

end

em1 = mean(err1);
ev1 = std(err1);
em2 = mean(err2);
ev2 = std(err2);

figure
hold on
errorbar(em1, ev1)
errorbar(em2, ev2)
%plot(1:Nout, model2.err(:,3), 'r-', 1:Nout, model1.err(:,3), 'b-');
legend('Batch HP', 'Stochastic HP')