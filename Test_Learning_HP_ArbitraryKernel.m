%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Comparison of different learning methods of Hawkes processes
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
addpath('Simulation')
addpath('Learning')



options.N = 200; % the number of sequences
options.Nmax = 100; % the maximum number of events per sequence
options.Tmax = 50; % the maximum size of time window
options.tstep = 0.1;
options.dt = 0.1;
options.M = 250;
options.GenerationNum = 10;
D = 3; % the dimension of Hawkes processes
nTest = 1;
nSeg = 5;
nNum = options.N/nSeg;


disp('Approximate simulation of Hawkes processes via branching process')
disp('Complicated gaussian kernel')
para1.kernel = 'gauss';
para1.w = 1.5; 
para1.landmark = 0:4:12;
L = length(para1.landmark);
para1.mu = rand(D,1)/D;
para1.A = zeros(D, D, L);
for l = 1:L
    para1.A(:,:,l) = (0.5^l)*(0.5+rand(D));
end
para1.A = 0.9*para1.A./max(abs(eig(sum(para1.A,3))));
para1.A = reshape(para1.A, [D, L, D]);
Seqs1 = Simulation_Branch_HP(para1, options);
%Seqs1 = Simulation_Thinning_HP(para1, options);


%%
disp('Learning Hawkes processes via different methods')
alg1.LowRank = 0;
alg1.Sparse = 1;
alg1.alphaS = 1;
alg1.GroupSparse = 0;
alg1.outer = 5;
alg1.rho = 0.1;
alg1.inner = 8;
alg1.thres = 1e-5;
alg1.Tmax = [];

alg2.alpha = 10000;
alg2.inner = 3;
alg2.inner_g = 100;
alg2.outer = 8;
alg2.thres = 1e-5;
alg2.Tmax = [];
Err = zeros(nTest, nSeg);

for n = 1:nTest
    for i = nSeg
       
        
        [A, Phi] = ImpactFunc( para1, options );
        
        disp('Maximum likelihood estimation and basis representation')        
        model1 = Initialization_Basis(Seqs1);
        model1 = Learning_MLE_Basis( Seqs1(1:i*nNum), model1, alg1 ); 
        [A1, Phi1] = ImpactFunc( model1, options );
        
        disp('Maximum likelihood estimation and ODE') 
        model2.M = 1000;
        model2.D = 2;
        model2.dt = 0.02;
        model2.g = rand(model2.M, model2.D);
        model2.g = model2.g./repmat(sum(model2.g),[model2.M,1]);
        model2.A = rand(D, model2.D, D)./(model2.D*D^2);
        model2.mu = rand(D,1)./D;
        model2 = Learning_MLE_ODE( Seqs1(1:i*nNum), model2, alg2 ); 
        [A2, Phi2] = ImpactFunc_ODE( model2 );
        
        disp('Least squares and discretization')
        model3.D = D;
        model3.h = 1;
        model3.k = floor(options.M * options.dt/model3.h);
        model3 = Initialization_Discrete(Seqs1(1:i*nNum));
        model3 = Learning_LS_Discrete( Seqs1(1:i*nNum), model3 );

        
        figure
        title('Complicated Gaussian kernels')
        for u = 1:D
            for v = 1:D
                subplot(D,D,D*(u-1)+v)
                hold on
                plot(options.dt*(0:(size(Phi,2)-1)), Phi(v,:,u), '-')
                plot(options.dt*(0:(size(Phi1,2)-1)), Phi1(v,:,u), '-')
                plot(options.dt*(0:(size(Phi2,2)-1)), Phi2(v,:,u), '-')
                plot(model3.h*(0:(size(model3.A,2)-1)), model3.A(v,:,u), '-')
                hold off
                axis tight
                legend('Real', 'MLE-Basis', 'MLE-ODE', 'LS')
                xlabel('Time interval between events')
                ylabel(['\phi', sprintf('%d%d', u, v)])
            end
        end
                
    end
end

