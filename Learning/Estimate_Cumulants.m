function model = Estimate_Cumulants(Seqs, model)

tic
Lambda = zeros(model.D,1);
C = zeros(model.D);
K = zeros(model.D);%, model.D, model.D);
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


for n = 1:length(Seqs)
    Ktmp = zeros(model.D);
    T = Seqs(n).Time(end)+eps;
    for i = 1:model.D
        indi = find(Seqs(n).Mark==i); 
        timei = Seqs(n).Time(indi);
        if ~isempty(indi)
            for m = 1:length(indi)
                indi1 = find(Seqs(n).Mark==i...
                    & Seqs(n).Time>=timei(m)-model.H ...
                    & Seqs(n).Time<=timei(m)+model.H);
                tmpi = length(indi1) - 2*model.H * Lambda(i);
            
                for j = 1:model.D
                    indj = find(Seqs(n).Mark==j);
                    timej = Seqs(n).Time(indj);
                    
                    dt = repmat(timei(:), [1, length(timej)]) - ...
                        repmat(timej(:)', [length(timei), 1]);
                    %dt = dt.*double(dt>0);
                    dt = 2*model.H - abs(dt);
                    dt = dt.*double(dt>0);
                    
                    indj1 = find(Seqs(n).Mark == j ...
                        & Seqs(n).Time>=timei(m)-model.H ...
                        & Seqs(n).Time<=timei(m)+model.H);
                    tmpj = length(indj1) - 2*model.H*Lambda(j);
                    
                    
                    Ktmp(i,j) = Ktmp(i,j) + tmpi * tmpj/T -...
                        Lambda(i)/T * sum(dt(:));
                end
            end
        end
    end
    K = K + Ktmp;
    if mod(n, 10)==0
        fprintf('Calculate 3rd order cumulant: n=%d/%d, time=%.2fsec\n',...
            n, length(Seqs), toc);
    end
end
K = K./length(Seqs) + 4 * model.H^2 * (Lambda(:).^2) * Lambda(:)'; 




% id = randperm(length(Seqs));
% for n = 1:min([length(id), model.Upper])
%     Ctmp = zeros(model.D);
%     Ktmp = zeros(model.D);
%     
%     T = Seqs(id(n)).Time(end)+eps;
%     for i = 1:model.D
%         indi = find(Seqs(id(n)).Mark==i);
%         if ~isempty(indi)
%             for m = 1:length(indi)
%                 indi1 = find(Seqs(id(n)).Mark==i ...
%                         & Seqs(id(n)).Time>=Seqs(id(n)).Time(indi(m))-model.H ...
%                         & Seqs(id(n)).Time<Seqs(id(n)).Time(indi(m))+model.H);
%                 tmpi = length(indi1) - 2*model.H * Lambda(i);
%                 
%                 for j = 1:model.D
%                     indj = find(Seqs(id(n)).Mark==j ...
%                         & Seqs(id(n)).Time>=Seqs(id(n)).Time(indi(m))-model.H ...
%                         & Seqs(id(n)).Time<Seqs(id(n)).Time(indi(m))+model.H);
%                     tmpj = length(indj) - 2*model.H*Lambda(j);
%                     Ctmp(i,j) = Ctmp(i,j) + tmpj;
%                     Ktmp(i,j) = Ktmp(i,j) + tmpi * tmpj;
%                 end
%             end
%         end
%     end
%     Ctmp = Ctmp./T;
%     Ktmp = Ktmp./T;
%     C = C + Ctmp;
%     K = K + Ktmp;
% end
% C = C/min([length(id), model.Upper]);
% 
% 
% for n = 1:min([length(id), model.Upper])
%     Ktmp = zeros(model.D);
%     T = Seqs(id(n)).Time(end)+eps;
%     
%     for j = 1:model.D
%         indj = find(Seqs(id(n)).Mark==j);
%         if ~isempty(indj)
%             for m = 1:length(indj)
%                 for i = 1:model.D
%                     indi1 = find(Seqs(id(n)).Mark==i ...
%                         & Seqs(id(n)).Time>=Seqs(id(n)).Time(indj(m))-2*model.H ...
%                         & Seqs(id(n)).Time<Seqs(id(n)).Time(indj(m))+2*model.H);
%                     indi2 = find(Seqs(id(n)).Mark==i ...
%                         & Seqs(id(n)).Time>=Seqs(id(n)).Time(indj(m))-2*model.H ...
%                         & Seqs(id(n)).Time<Seqs(id(n)).Time(indj(m)));
%                     tmpi1 = length(indi1) - 4*model.H*Lambda(i);
%                     Ktmp(i,j) = Ktmp(i,j) - 2*model.H*Lambda(i)/T * tmpi1;
% 
%                     if ~isempty(indi2)
%                         tmpi2 = sum(Seqs(id(n)).Time(indj(m)) - Seqs(id(n)).Time(indi2));
%                         Ktmp(i,j) = Ktmp(i,j) + 2*Lambda(i)/T*tmpi2;
%                     end
%                 end
%             end
%         end
%     end
%     
%     K = K + Ktmp;
% end
% 
% K = K/min([length(id), model.Upper]) - 4 * model.H^2 * (Lambda(:).^2) * Lambda(:)';

model.Cumulant1 = Lambda;
model.Cumulant2 = C;
model.Cumulant3 = K;