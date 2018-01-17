%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Superposed Hawkes Processes: various source_num, D = 10
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
%run setup
D = 10;
K = D;
N = 5;
nTest = 10;
rate = 0;

Err = zeros(2, 2, N, nTest);
ErrS = Err;
para = cell(K, 2, nTest);
Seqs = cell(K, 2, nTest);
models = cell(2, 2, N, nTest);


alg.LowRank = 0;
alg.Sparse = 0;
alg.GroupSparse = 0;
alg.outer = 10;
alg.rho = 0.1;
alg.inner = 10;
alg.thres = 2e-2;
alg.Tmax = [];

algmu = alg;
algmu.Sparse = 1;
algmu.alphaS = 0.1;

options.N = 300; % the number of sequences
options.Nmax = 100; % the maximum number of events per sequence
options.Tmax = 20; % the maximum size of time window
options.tstep = 0.1;
options.dt = 0.1;
options.M = 250;
options.GenerationNum = 200;


for n = 1:nTest
    A = rand(D, D);
    A(A<0.4) = 0;
    A = A + eye(D);
    A = 0.5*A./max(abs(eig(A)));
    AA = reshape(A, [D, 1, D]);
    L = zeros(K,2);
    for k = 1:K
        para{k, 1, n}.kernel = 'exp';
        para{k, 1, n}.landmark = 0;
        para{k, 1, n}.w = 1;
        para{k, 1, n}.mu = zeros(D,1);
        para{k, 1, n}.mu(k) = 0.5;
        para{k, 1, n}.A = reshape(A, [D, 1, D]);
        Seqs{k, 1, n} = SimulationFast_Thinning_ExpHP(para{k, 1, n}, options); 
        L(k,1) = length(Seqs{k, 1, n});
        
        para{k, 2, n}.kernel = 'exp';
        para{k, 2, n}.landmark = 0;
        para{k, 2, n}.w = 1;
        para{k, 2, n}.mu = ones(D,1)./(2*D);
        para{k, 2, n}.A = reshape(A, [D, 1, D]);
        Seqs{k, 2, n} = SimulationFast_Thinning_ExpHP(para{k, 2, n}, options); 
        L(k,2) = length(Seqs{k, 2, n});
    end
    
    
    model.A = rand(D,1,D)./(D^2);
    model.mu = rand(D,1)./D;
    model.kernel = 'exp';
    model.w = 1;
    model.landmark = 0;

        
        
    for nn = 1:N

        num = min([nn*(options.N/N); L(:)]);
        seq_sum1 = Seqs{1, 1, n}(1:num);
        seq_super1 = Seqs{1, 1, n}(1:num);
        seq_sum2 = Seqs{1, 2, n}(1:num);
        seq_super2 = Seqs{1, 2, n}(1:num);
        for kk = 2:K
            seq_sum1 = [seq_sum1, Seqs{kk, 1, n}(1:num)];
            seq_super1 = SuperPosition(seq_super1, ...
                                         Seqs{kk, 1, n}(1:num), rate);
            seq_sum2 = [seq_sum2, Seqs{kk, 2, n}(1:num)];
            seq_super2 = SuperPosition(seq_super2, ...
                                         Seqs{kk, 2, n}(1:num), rate);
        end            
         
        % multi-task Hawkes
        model1 = model;
%         model1.mu = 0.5*ones(1, num);
%         for k = 2:K
%             model1.mu = blkdiag(model1.mu, 0.5*ones(1, num));
%         end
        model1.mu = rand(D, num*K)./D;
        models{1, 1, nn, n} = model1;        
        models{1, 1, nn, n} = Learning_MLE_Basis_MTmu( ...
            seq_sum1, models{1, 1, nn, n}, algmu, 1 );

        model2 = model;
        %model2.mu = ones(D,num*K)./(2*D);%
        model2.mu = rand(D, num*K)./D;
        models{1, 2, nn, n} = model2;        
        models{1, 2, nn, n} = Learning_MLE_Basis_MTmu( ...
            seq_sum2, models{1, 2, nn, n}, algmu, 1 );
        
        
        % superposed Hawkes
        models{2, 1, nn, n} = model;
        %models{2, 1, nn, n}.mu = sum(model1.mu, 2)./num;
        models{2, 1, nn, n} = Learning_MLE_Basis2( ...
            seq_super1, models{2, 1, nn, n}, alg, 1 );
                
        models{2, 2, nn, n} = model;
        %models{2, 2, nn, n}.mu = sum(model2.mu, 2)./num;
        models{2, 2, nn, n} = Learning_MLE_Basis2( ...
            seq_super2, models{2, 2, nn, n}, alg, 1 );
        
            
        Err(1, 1, nn, n) = norm(models{1, 1, nn, n}.A(:) - ...
                AA(:))/norm(AA(:));
            
        Err(1, 2, nn, n) = norm(models{1, 2, nn, n}.A(:) - ...
                AA(:))/norm(AA(:));
            
        Err(2, 1, nn, n) = norm(models{2, 1, nn, n}.A(:) - ...
                AA(:))/norm(AA(:));
            
        Err(2, 2, nn, n) = norm(models{2, 2, nn, n}.A(:) - ...
                AA(:))/norm(AA(:));
            
            
        
    end
end

save('Result_superpose_error_syn_D10.mat', 'Err', 'models', ...    
                                        'para', 'options', 'K');