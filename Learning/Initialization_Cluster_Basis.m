function model = Initialization_Cluster_Basis(Seqs, ClusterNum, ...
                        baseType, bandwidth, landmark)

N = length(Seqs);
D = zeros(N,1);
for i = 1:N          
    D(i) = max(Seqs(i).Mark);
end
D = max(D);
model.K = ClusterNum;
model.D = D;

switch nargin
    case 2
        sigma = zeros(length(Seqs), 1);
        Tmax = zeros(length(Seqs), 1);

        for i = 1:length(Seqs)
            sigma(i) = ((4*std(Seqs(i).Time)^5)/(3*length(Seqs(i).Time)))^0.2; 
            Tmax(i) = Seqs(i).Time(end) + eps;
        end
        Tmax = mean(Tmax);
        
        model.kernel = 'gauss';
        model.w = mean(sigma);
        model.landmark = model.w*(0:ceil(Tmax/model.w));
        
    case 3
        model.kernel = baseType;
        model.w = 1;
        model.landmark = 0;
    case 4
        model.kernel = baseType;
        model.w = bandwidth;
        model.landmark = 0;
    otherwise
        model.kernel = baseType;
        model.w = bandwidth;
        model.landmark = landmark;
end


model.alpha = 1;
M = length(model.landmark);
model.beta = ones(D, M, model.K, D)./(M*D^2); 
% hyperparameter of Rayleigh prior of basic intensity (upper bound)
model.b = ones(D, model.K)./(D);  


% initialize label and responsibility randomly
label = ceil(model.K*rand(1,N));
model.R = full(sparse(1:N,label,1,N,model.K,N));
