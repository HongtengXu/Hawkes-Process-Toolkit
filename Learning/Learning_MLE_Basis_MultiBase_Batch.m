function model = Learning_MLE_Basis_MultiBase_Batch( Seqs, model, alg )
                                                        
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

if alg.storeLL
    model.LL = zeros(alg.epoch, 1);
end
if alg.storeErr
    model.err = zeros(alg.epoch+1, 3);
end

if alg.storeErr
    Err = zeros(1,3);
    Err(1) = norm(model.mu(:) - alg.truth.mu(:))/norm(alg.truth.mu(:));
    Err(2) = norm(model.A(:) - alg.truth.A(:))/norm(alg.truth.A(:));
    Err(3) = norm([model.mu(:); model.A(:)]-[alg.truth.mu(:); alg.truth.A(:)])...
        /norm([alg.truth.mu(:); alg.truth.A(:)]);
    model.err(1,:) = Err;
end

% initial 
if alg.LowRank
    UL = zeros(size(model.A));
    ZL = model.A;
end

if alg.Sparse
    US = zeros(size(model.A));
    ZS = model.A;
end

if alg.GroupSparse
    UG = zeros(size(model.A));
    ZG = model.A;
end

if alg.LowRankM
    UmuL = zeros(size(model.mu));
    ZmuL = model.mu;
end

if alg.SparseM
    UmuS = zeros(size(model.mu));
    ZmuS = model.mu;
end


D = size(model.A, 1);


tic;
for o = 1:alg.epoch
    
    rho = alg.rho * (1.1^o);
    
    for n = 1:alg.inner
        
        NLL = 0; % negative log-likelihood
        
        Amu = zeros(size(model.mu));
        Bmu = Amu;
        Cmu = Amu;        
        if alg.LowRankM
            Bmu = Bmu + rho*(UmuL-ZmuL);
            Amu = Amu + rho;
        end
        if alg.SparseM
            Bmu = Bmu + rho*(UmuS-ZmuS);
            Amu = Amu + rho;
        end        
        
        CmatA = zeros(size(model.A));
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
                        double(model.A(ui,:,:)>0).*repmat( GK(i,:), [1,1,D] );

                    ti = Time(i);             

                    lambdai = model.mu(ui,c);
                    pii = model.mu(ui,c);
                    pij = [];


                    if i>1

                        tj = Time(1:i-1);
                        uj = Event(1:i-1);

                        dt = ti - tj;
                        gij = Kernel(dt, model);
                        auiuj = model.A(uj, :, ui);
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

                NLL = NLL + (Tstop-Tstart).*sum(model.mu(:,c));
                NLL = NLL + sum( sum( GK.*sum(model.A(Event,:,:),3) ) );

            
            else
                warning('Sequence %d is empty!', c)
            end
        end
                
        % M-step: update parameters
        if alg.Sparse==0 && alg.GroupSparse==0 && alg.LowRank==0
            A = -CmatA./BmatA;            
        else            
            A = ( -BmatA + sqrt(BmatA.^2 - 4*AmatA.*CmatA) )./(2*AmatA);            
        end
        A(isnan(A))=0;
        A(isinf(A))=0;
            
        if alg.SparseM==0 && alg.LowRankM==0
            mu = Cmu./Bmu;%( -BA+sqrt(BA.^2-4*AA*CA) )./(2*AA);            
        else            
            mu = ( -Bmu + sqrt(Bmu.^2 + 4*Amu.*Cmu) )./(2*Amu);            
        end
        mu(isnan(mu))=0;
        mu(isinf(mu))=0;
        
        Error = norm([model.mu(:); model.A(:)]-[mu(:); A(:)])...
                        /norm([model.mu(:); model.A(:)]);
        model.A = A;
        model.mu = mu;
                    
    end
    fprintf('BatchOpt MHPs: Outer=%d, Obj=%f, RelErr=%f, Time=%0.2fsec\n',...
                o, NLL, Error, toc);
    
    % calculate Loglikelihood
    if alg.storeLL
        Loglike = Loglike_Basis( Seqs, model, alg );
        model.LL(o) = Loglike;
    end
    % calculate error
    if alg.storeErr
        Err = zeros(1,3);
        Err(1) = norm(model.mu(:) - alg.truth.mu(:))/norm(alg.truth.mu(:));
        Err(2) = norm(model.A(:) - alg.truth.A(:))/norm(alg.truth.A(:));
        Err(3) = norm([model.mu(:); model.A(:)]-[alg.truth.mu(:); alg.truth.A(:)])...
            /norm([alg.truth.mu(:); alg.truth.A(:)]);
        model.err(o+1,:) = Err;
    end
    
    
    if alg.LowRankM
        threshold = alg.betaLR/rho;
        ZmuL = SoftThreshold_LR( model.mu+UmuL, threshold );
        UmuL = UmuL + (model.mu-ZmuL);
    end

    if alg.SparseM
        threshold = alg.betaS/rho;
        ZmuS = SoftThreshold_S( model.mu+UmuS, threshold );
        UmuS = UmuS + (model.mu-ZmuS);
    end
    
    
    if alg.LowRank
        threshold = alg.alphaLR/rho;
        ZL = SoftThreshold_LR( model.A+UL, threshold );
        UL = UL + (model.A-ZL);
    end
    
    if alg.Sparse
        threshold = alg.alphaS/rho;
        ZS = SoftThreshold_S( model.A+US, threshold );
        US = US + (model.A-ZS);
    end

    if alg.GroupSparse
        threshold = alg.alphaGS/rho;
        ZG = SoftThreshold_GS( model.A+UG, threshold );
        UG = UG + (model.A-ZG);
    end    

                
end


