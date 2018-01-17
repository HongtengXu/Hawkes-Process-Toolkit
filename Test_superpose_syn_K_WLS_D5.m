%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Superposed Hawkes Processes: various source_num, D = 5
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
run setup
D = 5;
K = [2, 5, 10];
N = 5;
nTest = 10;
rate = 0;

Err = zeros(4, length(K), N, nTest);

para = cell(length(K), nTest);
Seqs = cell(K(end), nTest);
models = cell(4, length(K), N, nTest);

alg.LowRank = 0;
alg.Sparse = 0;
alg.GroupSparse = 0;
alg.outer = 10;
alg.rho = 0.1;
alg.inner = 10;
alg.thres = 1e-2;
alg.Tmax = [];
alg.w = 1;

algmu = alg;
algmu.Sparse = 1;
algmu.alphaS = 0.1;
algmu.w = 1;

options.N = 20; % the number of sequences
options.Nmax = 200; % the maximum number of events per sequence
options.Tmax = 100; % the maximum size of time window
options.tstep = 0.1;
options.dt = 0.1;
options.M = 250;
options.GenerationNum = 200;


for n = 1:nTest
    A = rand(D, D);
    A(A<0.4) = 0;
    A = A + eye(D);
    A = 0.5*A./max(abs(eig(A)));
    L = zeros(K(end),1);
    for k = 1:K(end)
        para{k, n}.kernel = 'exp';
        para{k, n}.landmark = 0;
        para{k, n}.w = 1;
        para{k, n}.mu = zeros(D,1);
        id = randperm(D);
        para{k, n}.mu(id(1)) = 0.1+rand;
        para{k, n}.A = reshape(A, [D, 1, D]);
        Seqs{k, n} = SimulationFast_Thinning_ExpHP(para{k, n}, options); 
        L(k) = length(Seqs{k, n});
    end
    
    for k = 1:length(K)%:-1:1
        model.A = rand(D,1,D)./(D^2);
        model.mu = rand(D,1)./D;
        model.kernel = 'exp';
        model.w = 1;
        model.landmark = 0;

        
        
        for nn = 1:N

            num = min([nn*(options.N/N); L]);
            seq_sum = Seqs{1, n}(1:num);
            seq_super = Seqs{1,n}(1:num);
            for kk = 2:K(k)
                seq_sum = [seq_sum, Seqs{kk, n}(1:num)];
                seq_super = SuperPosition(seq_super, ...
                                             Seqs{kk, n}(1:num), rate);
            end
            
%             seq_super = [];
%             for m = 1:K(k)
%                 seq_super_tmp = Seqs{1, n}(1:num);
%                 for kk = 2:K(k)
%                     tmp = Seqs{kk, n}(1:num);
%                     index = randperm(num);
%                     seq_super_tmp = SuperPosition(seq_super_tmp, ...
%                                             Seqs{kk, n}(index), rate);
%                 end
%                 seq_super = [seq_super, seq_super_tmp];
%             end
            
            models{1, k, nn, n} = model;        
            models{1, k, nn, n} = Learning_WLS_HP( ...
                Seqs{1, n}(1:num), models{1, k, nn, n}, alg, 'singlesource' );
            
            models{2, k, nn, n} = model;        
            models{2, k, nn, n} = Learning_WLS_HP( ...
                seq_sum, models{2, k, nn, n}, alg, 'singlesource' );
            
            model2 = model;
            model2.mu = rand(D, num*K(k))./D;
            models{3, k, nn, n} = model2;        
            models{3, k, nn, n} = Learning_WLS_HP( ...
                seq_sum, models{3, k, nn, n}, algmu, 'multisource' );
            
            models{4, k, nn, n} = model;  
            alg2 = alg;
            alg2.w = 1/sqrt(K(k));
            models{4, k, nn, n} = Learning_WLS_HP( ...
                seq_super, models{4, k, nn, n}, alg2, 'singlesource' );
                

            
            Err(1, k, nn, n) = norm(models{1, k, nn, n}.A(:) - ...
                A(:))/norm(A(:));
            
            Err(2, k, nn, n) = norm(models{2, k, nn, n}.A(:) - ...
                A(:))/norm(A(:));
            
            Err(3, k, nn, n) = norm(models{3, k, nn, n}.A(:) - ...
                A(:))/norm(A(:));
            
            Err(4, k, nn, n) = norm(models{4, k, nn, n}.A(:) - ...
                A(:))/norm(A(:));
            
%             A1 = models{1,k,nn,n}.A;
%             A2 = models{2,k,nn,n}.A;
%             imagesc([A,reshape(A1,[D,D]),reshape(A2,[D,D])]);
%             axis tight;axis equal;colorbar
            
        end
    end
end

save('Result_superpose_error_syn_WLS_D5.mat', 'Err', 'models', ...    
'para', 'options', 'K');