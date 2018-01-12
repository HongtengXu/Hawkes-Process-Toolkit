function model = Learning_MLE_Basis_Stoc( Seqs, model, alg )
                                                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Learning Hawkes processes via maximum likelihood estimation
% Different regularizers (low-rank, sparse, group sparse) of parameters and
% their combinations are considered, which are solved via StochasticADMM.
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
% KerFint = struct('dGK', []);
% KerF = struct('gij', []);

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

if alg.storeLL
    model.LL = zeros(alg.epoch, 1);
end
if alg.storeErr
    model.err = zeros(alg.epoch, 3);
end


tic;
for o = 1:alg.epoch % the number of epoch for stochastic opt
    
    rho = alg.rho * (1.1^o);
    
    idx_seq = randperm(length(Seqs));
    if length(Seqs)<alg.seqbatch;
        N_seq_batch = 1;
        Ns = length(Seqs);
    else
        N_seq_batch = floor(length(Seqs)/alg.seqbatch);
        Ns = alg.seqbatch;
    end
    
    for n = 1:N_seq_batch
        
        % initial para       
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
        for ns = 1:Ns
            c = idx_seq(Ns*(n-1)+ns); % the id of current sequence
            
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
                
                Nc = length(Time);
%                 if o == 1
%                     KerFint(c).dGK = zeros(alg.historyL, length(model.landmark), Nc-1);
%                     KerF(c).gij = zeros(alg.historyL, length(model.landmark), Nc-1);
%                 end
                
                idx_event = randperm(Nc);
                if Nc<alg.eventbatch
                    N_event_batch = 1;
                    Ne = Nc;
                else
                    N_event_batch = alg.eventbatch;
                    Ne = floor(Nc/alg.eventbatch);
                end
                TimeNew = [Tstart, Time];
                for m = 1:N_event_batch
                    for ne = 1:Ne
                        i = idx_event(Ne*(m-1)+ne);
                        ui = Event(i);
                        ti = Time(i);    
                        ti1 = TimeNew(i);
                        
                        
                        Amu = Amu + ti - ti1;
                        
                        
                        % calculate intensity
                        lambdai = muest(ui);
                        pii = muest(ui);
                        pij = [];

                        if i>1
                            start = max([1, i-alg.historyL]);
                            tj = Time(start:i-1);
                            uj = Event(start:i-1);

                            %if o==1 || isempty(KerF(c).gij)
                                dt = ti - tj;
                                dt1 = ti1 - tj;
                                gij = Kernel(dt, model);
                                
                                %KerF(c).gij(1:length(dt),:,i-1) = gij;
                                
                                GK = Kernel_Integration(dt, model);
                                GK1 = Kernel_Integration(dt1, model);
                                dGK = GK - GK1;
                                %KerFint(c).dGK(1:length(dt),:,i-1) = dGK;
                            %end
                            
                            auiuj = Aest(uj, :, ui);
                            %pij = auiuj .* KerF(c).gij(1:length(uj),:,i-1);%
                            pij = auiuj .* gij;
                            lambdai = lambdai + sum(pij(:));
                        end
                        pii = pii./lambdai;
                        
                        
                        if i>1
                            pij = pij./lambdai;
                            if ~isempty(pij) && sum(pij(:))>0
                                for j = 1:length(uj)
                                    uuj = uj(j);
                                    CmatA(uuj,:,ui) = CmatA(uuj,:,ui) - pij(j,:);
                                    
                                    BmatA(uuj,:,:) = BmatA(uuj,:,:) + ...
                                        double(Aest(ui,:,:)>0).*...
                                        repmat( dGK(j,:), [1,1,D] );
                                        %repmat( KerFint(c).dGK(j,:,i-1), [1,1,D] );
                                end
                            end
                        end
                        Bmu(ui) = Bmu(ui) + pii;
                        
                    end
                end
                
%                 Amu = Amu + Tstop - Tstart;
% 
%                 dT = Tstop - Time;
%                 GK = Kernel_Integration(dT, model);
% 
%                 Nc = length(Time);
%                 idx_event = randperm(Nc);
% 
%                 if Nc<alg.eventbatch
%                     N_event_batch = 1;
%                     Ne = Nc;
%                 else
%                     N_event_batch = alg.eventbatch;
%                     Ne = floor(Nc/alg.eventbatch);
%                 end
%                 
%                 for m = 1:N_event_batch
%                     for ne = 1:Ne
%                         i = idx_event(Ne*(m-1) + ne);
%                         ui = Event(i);
% 
%                         BmatA(ui,:,:) = BmatA(ui,:,:)+...
%                             double(Aest(ui,:,:)>0).*repmat( GK(i,:), [1,1,D] );
% 
%                         ti = Time(i);             
%                         lambdai = muest(ui);
%                         pii = muest(ui);
%                         pij = [];
% 
%                         if i>1
%                             start = max([1, i-alg.historyL]);
%                             tj = Time(start:i-1);
%                             uj = Event(start:i-1);
% 
%                             dt = ti - tj;
%                             gij = Kernel(dt, model);
%                             auiuj = Aest(uj, :, ui);
%                             pij = auiuj .* gij;
%                             lambdai = lambdai + sum(pij(:));
%                         end
% 
%                         pii = pii./lambdai;
% 
%                         if i>1
%                             pij = pij./lambdai;
%                             if ~isempty(pij) && sum(pij(:))>0
%                                 for j = 1:length(uj)
%                                     uuj = uj(j);
%                                     CmatA(uuj,:,ui) = CmatA(uuj,:,ui) - pij(j,:);
%                                 end
%                             end
%                         end
% 
%                         Bmu(ui) = Bmu(ui) + pii;
%                     end
%                 end           
            else
                warning('Sequence %d is empty!', c)
            end
        end
                
        % M-step: update parameters
        mu = Bmu./Amu;        
        if alg.Sparse==0 && alg.GroupSparse==0 && alg.LowRank==0
            A = -CmatA./BmatA;%( -BA+sqrt(BA.^2-4*AA*CA) )./(2*AA);
            A(isnan(A))=0;
            A(isinf(A))=0;
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
        fprintf('Epoch=%d, EventBatch=%d, RelErr=%f, Time=%0.2fsec\n',...
                o, n, Err, toc);
            
        if Err<alg.thres || (o==alg.epoch && n==N_seq_batch)
            break;
        end    
    end
    
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
        model.err(o,:) = Err;
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


