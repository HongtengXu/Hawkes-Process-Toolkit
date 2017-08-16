%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compare different clustering method of event sequences
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear
addpath('Simulation')
addpath('Learning')



options.N = 100; % the number of sequences
options.Nmax = 100; % the maximum number of events per sequence
options.Tmax = 50; % the maximum size of time window
options.tstep = 0.1;
options.dt = 0.1;
options.M = 250;
options.GenerationNum = 10;
D = 3; % the dimension of Hawkes processes
K = 2; % the number of clusters
nTest = 5;
nSeg = 5;
nNum = options.N/nSeg;


disp('Approximate simulation of Hawkes processes via branching process')
disp('Simple exponential kernel')
para1.kernel = 'exp';
para1.landmark = 0;
L = length(para1.landmark);
para1.mu = rand(D,1)/D;
para1.A = zeros(D, D, L);
for l = 1:L
    para1.A(:,:,l) = (0.7^l)*rand(D);
end
para1.A = 0.5 * para1.A./max(abs(eig(para1.A)));
para1.w = 0.5; 
Seqs1 = Simulation_Branch_HP(para1, options);

disp('Complicated gaussian kernel')
para2.kernel = 'gauss';
para2.landmark = 0:3:12;
L = length(para2.landmark);
para2.mu = rand(D,1)/D;
para2.A = zeros(D, D, L);
for l = 1:L
    para2.A(:,:,l) = (0.9^l)*rand(D);
end
para2.A = 0.25 * para2.A./max(abs(eig(sum(para2.A,3))));
para2.A = reshape(para2.A, [D, L, D]);
para2.w = 1; 
Seqs2 = Simulation_Branch_HP(para2, options);

SeqsMix = [Seqs1,Seqs2];

% ground truth
GT = blkdiag(ones(options.N), ones(options.N));

% Mixture of Hawkes processes
% initialize
model = Initialization_Cluster_Basis(SeqsMix, 2);
alg.outer = 8;
alg.rho = 0.1;
alg.inner = 5;
alg.thres = 1e-5;
alg.Tmax = [];
model = Learning_Cluster_Basis(SeqsMix, model, alg);
[~, labels1] = max(model.R, [],2);
Est1 = zeros(length(labels1));
for i = 1:length(labels1)
    for j = 1:length(labels1)
        if labels1(i)==labels1(j)
            Est1(i,j) = 1;
        end
    end
end

% Distance metric of marked point processes
configure.M = [options.Tmax; D];
configure.id = {1,2,[1,2]};
configure.tau = 1;
configure.W = 5;
configure.epoch = 15;
configure.lr = 1e-4;
[configure, obj] = Estimate_Weight(configure, SeqsMix);
Dis = zeros(length(SeqsMix));
tic
for n = 1:length(SeqsMix)-1
    Xn = [SeqsMix(n).Time; SeqsMix(n).Mark];
    for m = n+1:length(SeqsMix)
        Xm = [SeqsMix(m).Time; SeqsMix(m).Mark];
        [Dis(n,m), ~] = DistanceSum_MPP(Xn, Xm, configure);
    end
    toc
end
Dis = Dis + Dis';
SimilarMat = exp(-Dis.^2./(2*var(Dis(:))));


figure
subplot(131)
imshow(GT,[])
title('Ground truth of clusters')
subplot(132)
imshow(Est1,[])
title('Mixture of Hawkes processess')
subplot(133)
imshow(SimilarMat,[])
title('Distance metrics of MPP')
