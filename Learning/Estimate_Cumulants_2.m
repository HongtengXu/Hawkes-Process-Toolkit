function model = Estimate_Cumulants_2(Seqs, model)

tic
Lambda = zeros(model.D,1);
C = zeros(model.D);
for n = 1:length(Seqs)
    T = Seqs(n).Time(end)+eps;
    for i = 1:model.D
        ind = find(Seqs(n).Mark == i);
        Lambda(i) = Lambda(i) + length(ind)/T;
    end
end
Lambda = Lambda/length(Seqs);

for n = 1:length(Seqs)
    Ctmp = zeros(model.D);
    T = Seqs(n).Time(end)+eps;
    for i = 1:model.D
        indi = find(Seqs(n).Mark==i);
        timei = Seqs(n).Time(indi);
        if ~isempty(indi)
            for m = 1:length(indi)
                for j = 1:model.D
                    indj = find(Seqs(n).Mark == j ...
                        & Seqs(n).Time>=timei(m)-model.H ...
                        & Seqs(n).Time<=timei(m)+model.H);
                    Ctmp(i,j) = Ctmp(i,j) + length(indj) - 2*model.H*Lambda(j);
                end
            end
        end
    end
    C = C + Ctmp./T;
    if mod(n, 10)==0
        fprintf('Calculate 2nd order cumulant: n=%d/%d, time=%.2fsec\n',...
            n, length(Seqs), toc);
    end
end
C = C./length(Seqs);

model.Cumulant1 = Lambda;
model.Cumulant2 = C;
