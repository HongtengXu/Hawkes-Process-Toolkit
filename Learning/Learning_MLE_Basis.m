function model = Learning_MLE_Basis( Seqs, model, alg )
                                                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Learning Hawkes processes via maximum likelihood estimation
% Different regularizers (low-rank, sparse, group sparse) of parameters and
% their combinations are considered, which are solved via ADMM.
%
% Reference:
% Xu, Hongteng, Mehrdad Farajtabar, and Hongyuan Zha. 
% "Learning Granger Causality for Hawkes Processes." 
% International Conference on Machine Learning (ICML). 2016.
%
% Provider:
% Hongteng Xu @ Georgia Tech
% June. 10, 2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initial 
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


tic;
for o = 1:alg.outer
    
    rho = alg.rho * (1.1^o);
    
    for n = 1:alg.inner
        
        NLL = 0; % negative log-likelihood
        
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
        
        % E-step: evaluate the responsibility using the current parameters    
        for c = 1:length(Seqs)
            Time = Seqs(c).Time;
            Event = Seqs(c).Mark;
            if isempty(alg.Tmax)
                Tend = Time(end)+eps;
            else
                Tend = alg.Tmax;
                indt = Time < alg.Tmax;
                Time = Time(indt);
                Event = Event(indt);
            end

            Amu = Amu + Tend;
            
            dT = Tend - Time;
            GK = Kernel_Integration(dT, model);
            
            Nc = length(Time);
            
            for i = 1:Nc
                
                ui = Event(i);
                
                BmatA(ui,:,:) = BmatA(ui,:,:)+...
                    double(Aest(ui,:,:)>0).*repmat( GK(i,:), [1,1,D] );
                
                ti = Time(i);             
                    
                lambdai = muest(ui);
                pii = muest(ui);
                pij = [];
                          
                    
                if i>1
                    
                    tj = Time(1:i-1);
                    uj = Event(1:i-1);

                    dt = ti - tj;
                    gij = Kernel(dt, model);
                    auiuj = Aest(uj, :, ui);
                    pij = auiuj .* gij;
                    lambdai = lambdai + sum(pij(:));
                end

                NLL = NLL - log(lambdai);
                pii = pii./lambdai;
                
                if i>1
                    pij = pij./lambdai;
                    if ~isempty(pij) && sum(pij(:))>0
                        for j = 1:length(uj)
                            uuj = uj(j);
                            CmatA(uuj,:,ui) = CmatA(uuj,:,ui) - pij(j,:);
                        end
                    end
                end
                
                Bmu(ui) = Bmu(ui) + pii;
                
            end
            
            NLL = NLL + Tend.*sum(muest);
            NLL = NLL + sum( sum( GK.*sum(Aest(Event,:,:),3) ) );

            
            
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
                o, n, NLL, Err, toc);
            
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


