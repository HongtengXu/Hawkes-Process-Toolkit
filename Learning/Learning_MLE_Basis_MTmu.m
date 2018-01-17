function model = Learning_MLE_Basis_MTmu( Seqs, model, alg )
                                                        
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
    UL = zeros(size(muest));
    ZL = muest;
end

if alg.Sparse
    US = zeros(size(muest));
    ZS = muest;
end


D = size(Aest, 1);


tic;
for o = 1:alg.outer
    
    rho = alg.rho * (1.1^o);
    
    for n = 1:alg.inner
        
        NLL = 0; % negative log-likelihood
        
        Amu = zeros(size(muest));
        Bmu = Amu;
        Cmu = Amu;
        if alg.updatemu==1
            if alg.LowRank
                Bmu = Bmu + rho*(UL-ZL);
                Amu = Amu + rho;
            end
            if alg.Sparse
                Bmu = Bmu + rho*(US-ZS);
                Amu = Amu + rho;
            end
        end
        
        
        CmatA = zeros(size(Aest));
        AmatA = CmatA;
        BmatA = CmatA;
        
        
        % E-step: evaluate the responsibility using the current parameters    
        for c = 1:length(Seqs)
            if ~isempty(Seqs(c).Time)
                Time = Seqs(c).Time;
                Event = Seqs(c).Mark;
                Tstart = Seqs(c).Start;

                if isempty(alg.Tmax)
                    Tstop = Seqs(c).Stop;
                else
                    Tstop = alg.Tmax;
                    indt = Time < alg.Tmax;
                    Time = Time(indt);
                    Event = Event(indt);
                end

                
                Bmu(:,c) = Bmu(:,c) + Tstop - Tstart;

                dT = Tstop - Time;
                GK = Kernel_Integration(dT, model);

                Nc = length(Time);

                for i = 1:Nc

                    ui = Event(i);

                    BmatA(ui,:,:) = BmatA(ui,:,:)+...
                        double(Aest(ui,:,:)>0).*repmat( GK(i,:), [1,1,D] );

                    ti = Time(i);             

                    lambdai = muest(ui,c);
                    pii = muest(ui,c);
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

                    Cmu(ui, c) = Cmu(ui, c) + pii;

                end

                NLL = NLL + (Tstop-Tstart).*sum(muest(:,c));
                NLL = NLL + sum( sum( GK.*sum(Aest(Event,:,:),3) ) );

            
            else
                warning('Sequence %d is empty!', c)
            end
        end
                
        % M-step: update parameters
        
        A = -CmatA./BmatA;%( -BA+sqrt(BA.^2-4*AA*CA) )./(2*AA);
        A(isnan(A))=0;
        A(isinf(A))=0;
            
        if alg.updatemu==1
            if alg.Sparse==0 && alg.LowRank==0
                mu = Cmu./Bmu;%( -BA+sqrt(BA.^2-4*AA*CA) )./(2*AA);
                mu(isnan(mu))=0;
                mu(isinf(mu))=0;
            else            
                mu = ( -Bmu + sqrt(Bmu.^2 + 4*Amu.*Cmu) )./(2*Amu);
                mu(isnan(mu))=0;
                mu(isinf(mu))=0;
            end
        
        else
            mu = muest;
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
    
    if alg.updatemu==1
        if alg.LowRank
            threshold = alg.alphaLR/rho;
            ZL = SoftThreshold_LR( muest+UL, threshold );
            UL = UL + (muest-ZL);
        end

        if alg.Sparse
            threshold = alg.alphaS/rho;
            ZS = SoftThreshold_S( muest+US, threshold );
            US = US + (muest-ZS);
        end
    end
    
        

                
end


