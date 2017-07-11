function [Prob, Delta] = ChosenProbability(X, configure)

L = size(X,2);
Prob = zeros(L-configure.W+1);
Dis = Prob;
dis = zeros(length(configure.weight), L-configure.W+1, L-configure.W+1);
Delta = dis;

for t1 = 1:L-configure.W+1
    for t2 = 1:L-configure.W+1
        
        if t1~=t2
            X1 = X(:, t1:t1+configure.W-1);
            X2 = X(:, t2:t2+configure.W-1);
            [Dis(t1, t2), dis(:, t1, t2)] = DistanceSum_MPP(X1, X2, configure);
            Prob(t1, t2) = exp(-Dis(t1, t2)^2);
        end
    end
end


Prob = Prob./repmat(sum(Prob,2)+eps, [1, size(Prob,2)]);


for t1 = 1:L-configure.W+1
    for t2 = 1:L-configure.W+1        
        if t1~=t2
            Delta(:, t1, t2) = sum(repmat(Prob(t1, :) .* Dis(t1, :), ...
                [size(dis,1),1,1]) .* dis(:, :, t1), 2) ...
                - (size(Prob,1)-1).* Dis(t1, t2).* dis(:, t1, t2);            
        end
    end
end