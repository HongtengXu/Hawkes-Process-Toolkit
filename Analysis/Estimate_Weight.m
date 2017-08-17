function [configure, obj] = Estimate_Weight(configure, Seqs)


% initialization
configure.weight = rand(length(configure.id),1);
tau = configure.tau;
obj = zeros(configure.epoch * length(Seqs), 1);

tic
for n = 1:configure.epoch
    ind = randperm(length(Seqs));
    lr = configure.lr * (0.9)^(n-1);
    for m = 1:length(Seqs)
        X = [Seqs(ind(m)).Time; Seqs(ind(m)).Mark];
        [Prob, Delta] = ChosenProbability(X, configure);
        
        grad = 0;
        for t1 = 1:size(Prob,1)-tau
            for t2 = 1:size(Prob,2)-tau
                if t1~=t2
                    obj((n-1)*length(Seqs)+m) = obj((n-1)*length(Seqs)+m) +...
                       Prob(t1, t2)*Prob(t1+tau, t2+tau); 
                    grad = grad + ...
                        2*Prob(t1, t2)*Prob(t1+tau, t2+tau)*...
                        (Delta(:,t1,t2)+Delta(:,t1+tau,t2+tau));
                end
            end
        end        
        configure.weight = configure.weight - lr * grad;
        
        fprintf('epoch=%d, #seq=%d/%d, obj=%f, ||grad||=%.4f, time=%.2fsec\n',...
            n, m, length(Seqs), obj((n-1)*length(Seqs)+m), norm(grad), toc);
    end
        
end