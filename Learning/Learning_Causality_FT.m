function model = Learning_Causality_FT(Seqs, model)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Learning Granger causality of Hawkes processes via a time series-based
% method.
%
% Reference:
% Etesami, Jalal, et al. 
% "Learning network of multivariate Hawkes processes: 
% a time series approach." UAI, 2016.
%
% Provider:
% Hongteng Xu @ Georgia Tech
% June. 12, 2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = length(Seqs);
Lambda = zeros(model.D, 1);
Sigma = zeros(model.D, model.D, model.Tau);

for i = 1:length(Seqs)
    Time = Seqs(i).Time;
    Event = Seqs(i).Mark;
    Y = zeros(model.D, 1);
    for d = 1:model.D
        ind = find(Event == d);
        Y(d) = length(ind)/(Time(end)+eps);
    end
    Lambda = Lambda+(1/N)*Y;
end


for i = 1:length(Seqs)
    Time = Seqs(i).Time;
    Event = Seqs(i).Mark;
    Th = floor(Time(end)./model.h);
    X = zeros(model.D, Th);
    for d = 1:model.D
        ind = Event == d;
        Time_d = Time(ind);
        for t = 1:Th
            td = find(Time_d<=t*model.h);
            X(d, t) = length(td) - Lambda(d)*(t*model.h);
        end
    end
    
    dX1 = X(:,2:end) - X(:,1:end-1);
    for tau = 1:model.Tau
        dXt = X(:,1+tau:Th) - X(:,1:Th-tau);
        tmp = (1/(Time(end)+eps)) * dX1(:,1:Th-tau) * dXt';
        Sigma(:,:,tau) = 1/length(Seqs) * tmp;
    end
        
end


LOmega = zeros(model.D, model.D, model.N);
Spec = zeros(1, model.N);
for d1 = 1:model.D
    for d2 = 1:model.D
        sigma = Sigma(d1,d2,:);
        tmp = fft(sigma(:)', model.N);
        LOmega(d1, d2, :) = reshape(tmp, [1,1,model.N]);
        if d1==d2            
            Spec = Spec + 1./abs(tmp);
        end
    end
end
Spec = 1./Spec;
ind = find(Spec(1:round(end/2))==0);
if isempty(ind)
    [~, ind] = min(Spec);
end
beta = 2*pi * ind./model.N;
LOmega = LOmega(:,:,ind);

bL = zeros(model.D * length(ind), model.D);
AL = zeros(model.D * length(ind));

for i = 1:length(ind)
    bL(1+(i-1)*model.D:i*model.D,:) = LOmega(:,:,i);
    for j = 1:length(ind)
        AL(1+(i-1)*model.D:i*model.D, 1+(j-1)*model.D:i*model.D) = ...
            (diag(Lambda) + LOmega(:,:,j) + LOmega(:,:,i)')./(beta(i)+beta(j));
    end
end

A = abs(AL\bL);
model.GCG = zeros(model.D);
for i = 1:length(ind)
    model.GCG = model.GCG + A(1+(i-1)*model.D:i*model.D,:);
end
model.GCG = model.GCG;
