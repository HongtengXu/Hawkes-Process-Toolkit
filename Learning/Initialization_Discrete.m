function [model, Tmax] = Initialization_Discrete(Seqs, step, order)

switch nargin
    case 1
        Num = zeros(length(Seqs), 1);
        Tmax = zeros(length(Seqs), 1);
        for i = 1:length(Seqs)
            Num(i) = length(Seqs(i).Time);
            Tmax(i) = Seqs(i).Time(end) + eps;
        end
        model.h = mean(Tmax./Num);   
        model.k = floor(min(0.7*Tmax/model.h));
    case 2
        model.h = step;
        Tmax = zeros(length(Seqs), 1);
        for i = 1:length(Seqs)
            Tmax(i) = Seqs(i).Time(end) + eps;
        end
        model.k = floor(min(0.7*Tmax/model.h));
    case 3
        model.h = step;
        model.k = order;
end

D = zeros(length(Seqs),1);
for i = 1:length(Seqs)            
    D(i) = max(Seqs(i).Mark);
end
model.D = max(D);