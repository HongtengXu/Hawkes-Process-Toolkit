%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compare different Granger causality learning method
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear
addpath('Simulation')
addpath('Learning')



options.N = 2000; % the number of sequences
options.Nmax = 500; % the maximum number of events per sequence
options.Tmax = 500; % the maximum size of time window
options.tstep = 0.1;
options.dt = 0.1;
options.M = 250;
options.GenerationNum = 100;
D = 10; % the dimension of Hawkes processes



disp('Approximate simulation of Hawkes processes via branching process')
disp('Complicated gaussian kernel')
para1.kernel = 'exp';
para1.w = 2; 
para1.landmark = 0;%:4:12;
L = length(para1.landmark);
para1.mu = rand(D,1)/D;
para1.A = zeros(D, D, L);
for l = 1:L
    para1.A(:,:,l) = (0.5^l)*(0.5+ones(D));
end
%mask = double(triu(rand(D))>0);%
mask = rand(D).*double(rand(D)>0.7);
% mask = zeros(D);
% maks(1:D/2,1:D/2) = 1;
para1.A = para1.A.*repmat(mask, [1,1,L]);
para1.A = 0.25*para1.A./max(abs(eig(sum(para1.A,3))));
tmp = para1.A;
para1.A = reshape(para1.A, [D, L, D]);
for di = 1:D
    for dj = 1:D
        phi = tmp(di, dj, :);
        para1.A(dj, :, di) = phi(:);
    end
end
Seqs1 = Simulation_Branch_HP(para1, options);
%Seqs1 = Simulation_Thinning_HP(para1, options);


%%
[A, Phi] = ImpactFunc( para1, options );

disp('Learning Hawkes processes via different methods')



% disp('Maximum likelihood estimation and basis representation') 
% alg1.LowRank = 0;
% alg1.Sparse = 1;
% alg1.alphaS = 1;
% alg1.GroupSparse = 1;
% alg1.alphaGS = 100;
% alg1.outer = 5;
% alg1.rho = 0.1;
% alg1.inner = 8;
% alg1.thres = 1e-5;
% alg1.Tmax = [];
% model1 = Initialization_Basis(Seqs1);
% model1 = Learning_MLE_Basis( Seqs1, model1, alg1 ); 
% [A1, Phi1] = ImpactFunc( model1, options );

disp('Learning causality via cumulants of Hawkes processes')
model2.D = D;
model2.H = 1.5;%0;%options.Tmax/2;
model2.Upper = 1500;
model2.Rreal = inv(eye(D)-A);
model2.mureal = para1.mu;
alg2.iterNum = 100;
alg2.lr = 1e-3;
alg2.alpha = 0.5;
model2 = Learning_Causality_Cumulants(Seqs1, model2, alg2);

% disp('Learning causality via least squares')
% model3.D = D;
% model3.h = 1;
% model3.k = floor(options.M * options.dt/model3.h);
% model3 = Learning_LS_Discrete( Seqs1, model3 );
% % model3.D = D;
% % model3.Tau = 1;
% % model3.h = 1;
% % model3.N = 1024;
% % model3 = Learning_Causality_FT(Seqs1, model3);


figure
subplot(221)        
imagesc(A)
title('Ground truth of infectivity')
axis square
colorbar
subplot(223)        
imagesc(A>0)
%title('Estimated infectivity-MLE')
colorbar
axis square
subplot(222)        
imagesc(model2.GCG)
title('Estimated infectivity-Cumulant')
colorbar
axis square
subplot(224)        
imagesc(model2.GCG>0)
%title('Estimated infectivity-LS')
colorbar
axis square
        


