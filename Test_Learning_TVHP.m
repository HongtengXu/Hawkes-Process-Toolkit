%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Comparison of different learning methods of Hawkes processes
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
addpath('Simulation')
addpath('Learning')



options.N = 500; % the number of sequences
options.Nmax = 100; % the maximum number of events per sequence
options.Tmax = 100; % the maximum size of time window
options.type = 'Est';
options.M = 500;
options.dt = options.Tmax/options.M;

D = 3; % the dimension of Hawkes processes
nTest = 1;
nSeg = 5;
nNum = options.N/nSeg;


disp('Approximate simulation of Time-varying Hawkes processes')
disp('Complicated gaussian kernel')
para1.kernel = 'gauss';
para1.w = 3;
para1.wt = 1;
para1.landmark = 0:5:options.Tmax;
L = length(para1.landmark);
para1.mu = rand(D,1)/D;
para1.A = zeros(D, D, L);
for l = 1:L
    para1.A(:,:,l) = rand(D)./D;
end
%para1.A = 0.9*para1.A./max(abs(eig(sum(para1.A,3))));
para1.A = reshape(para1.A, [D, L, D]);
Seqs1 = Simulation_Thinning_TVHP(para1, options);


%%
disp('Learning Hawkes processes via different methods')
alg1.LowRank = 0;
alg1.Sparse = 1;
alg1.alphaS = 0.01;
alg1.GroupSparse = 0;
alg1.outer = 5;
alg1.rho = 0.1;
alg1.inner = 8;
alg1.thres = 1e-5;
alg1.Tmax = [];

Err = zeros(nTest, nSeg);

for n = 1:nTest
    for i = nSeg
       
        
        [A, Phi] = ImpactFunc( para1, options );
        
        disp('Maximum likelihood estimation and basis representation')        
        model1.kernel = 'gauss';
        model1.sigmaA = 3;
        model1.sigmaT = 1;
        model1.w = model1.sigmaT;
        model1.landmark = 0:5:options.Tmax;
        model1.A = rand(D, length(model1.landmark), D)./(D^2 * length(model1.landmark));
        model1.mu = rand(D, 1)./D;
        
        model1 = Learning_MLE_TVHP_Basis( Seqs1(1:i*nNum), model1, alg1 ); 
        [A1, Phi1] = ImpactFunc( model1, options );
        
        
        figure
        title('Complicated Gaussian kernels')
        for u = 1:D
            for v = 1:D
                subplot(D,D,D*(u-1)+v)
                hold on
                plot(options.dt*(0:(size(Phi,2)-1)), Phi(v,:,u), '-')
                plot(options.dt*(0:(size(Phi1,2)-1)), Phi1(v,:,u), '-')
                hold off
                axis tight
                legend('Real', 'MLE-Basis')
                xlabel('Time interval between events')
                ylabel(['\phi', sprintf('%d%d', u, v)])
            end
        end
                
    end
end

