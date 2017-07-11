function model = Learning_MLE_TVHP_Basis( Seqs, model, alg )

tic
% initialize
Aest = model.A;        
muest = model.mu;

if alg.LowRank
    UL = zeros(size(Aest));
    ZL = Aest;
end

if alg.Sparse
    US = zeros(size(Aest));
    ZS = Aest;
end

if alg.GroupSparse
    UG = zeros(size(Aest));
    ZG = Aest;
end

D = size(Aest, 1);


for o = 1:alg.outer
    rho = alg.rho * 1.1^o;
    
    % update model.A, mu
    for l = 1:alg.inner
        
        NLLinner = 0;
        Amu = zeros(D, 1);
        Bmu = Amu;
        
        
        CmatA = zeros(size(Aest));
        AmatA = CmatA;
        BmatA = CmatA;
        if alg.LowRank
            BmatA = BmatA + rho*(UL-ZL);
            AmatA = AmatA + rho;
        end
        if alg.Sparse
            BmatA = BmatA + rho*(US-ZS);
            AmatA = AmatA + rho;
        end
        if alg.GroupSparse
            BmatA = BmatA + rho*(UG-ZG);
            AmatA = AmatA + rho;
        end
        
        
        Amc = sum(model.A,3);
        for c = 1:length(Seqs)
            
            Time = Seqs(c).Time;
            Event = Seqs(c).Mark;
            if isempty(alg.Tmax)
                Tstop = Time(end)+eps;
                Tstart = max([Time(1)-eps, 0]);
            else
                Tstop = alg.Tmax;
                Tstart = max([Time(1)-eps, 0]);
                indt = Time < alg.Tmax;
                Time = Time(indt);
                Event = Event(indt);
            end

            
            Amu = Amu + (Tstop-Tstart);
            
            In = length(Times);
            kappat = Kernel(Time, model);
            
            
            for i = 1:In
                ti = Time(i);
                ci = Event(i);
                lambdai = muest(ci);
                pii = muest(ci);
                
                
                lower = ti;
                if i==In
                    upper = Tstop;
                else
                    upper = Time(i+1);
                end
                KGt = IntKernelComp( Time(1:i), model, upper, lower);
                for j = 1:i                        
                    BmatA(Event(j), :, :) = BmatA(Event(j), :, :) +...
                        repmat(KGt(:,j)', [1,1,D]);
                end
                NLLinner = NLLinner + sum(sum(KGt'.*Amc(Event(1:i),:)));
                
                if i>1
                    tj = Time(1:i-1);
                    cj = Event(1:i-1);
                    
                    acicj = model.A(cj,:,ci);
                    tmp.landmark = tj;
                    tmp.kernel = model.kernel;
                    tmp.w = model.w;
                    gt = Kernel(ti, tmp);               
                    
                    pijm = repmat(gt(:), [1, size(acicj,2)])...
                        .* acicj...
                        .* repmat(kappat(i,:), [i-1, 1]);
                    
                    lambdai = lambdai + sum(pijm(:));
                
                    pijm = pijm./lambdai;
                    
                    for j = 1:i-1
                        CmatA(cj(j), :, ci) = CmatA(cj(j), :, ci) -...
                            pijm(j, :);                        
                    end
                end
                
                NLLinner = NLLinner - log(lambdai);
                
                pii = pii./lambdai;
                Bmu(ci) = Bmu(ci)+pii;
                
                
            end
            
            NLLinner = NLLinner + (Tstop-Tstart)*sum(model.mu(:));
            
        end
        
        
        % M-step: update parameters
        mu = Bmu./Amu;        
        if alg.Sparse==0 && alg.GroupSparse==0 && alg.LowRank==0
            A = -CmatA./BmatA;%( -BA+sqrt(BA.^2-4*AA*CA) )./(2*AA);
        else            
            A = ( -BmatA + sqrt(BmatA.^2 - 4*AmatA.*CmatA) )./(2*AmatA);
            A(isnan(A))=0;
            A(isinf(A))=0;
        end
        
        
        
        % check convergence
        Err=sum(sum(sum(abs(A-Aest))))/sum(abs(Aest(:)));
        Aest = A;
        muest = mu;
        model.A = Aest;
        model.mu = muest;
        fprintf('Outer=%d, Inner=%d, Obj=%f, RelErr=%f, Time=%0.2fsec\n',...
                o, n, NLLinner, Err, toc);
            
        if Err<alg.thres || (o==alg.outer && n==alg.inner)
            break;
        end    
    end
    
    if alg.LowRank
        threshold = alg.alphaLR/rho;
        ZL = SoftThreshold_LR( Aest+UL, threshold );
        UL = UL + (Aest-ZL);
    end
    
    if alg.Sparse
        threshold = alg.alphaS/rho;
        ZS = SoftThreshold_S( Aest+US, threshold );
        US = US + (Aest-ZS);
    end

    if alg.GroupSparse
        threshold = alg.alphaGS/rho;
        ZG = SoftThreshold_GS( Aest+UG, threshold );
        UG = UG + (Aest-ZG);
    end
        

                
end


