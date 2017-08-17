function model = Learning_MLE_ODE( Seqs, model, alg )
                                                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Learning Hawkes processes via maximum likelihood estimation with the help
% of ordinary differential equations (ODE)
%
% Reference:
% Zhou, Ke, Hongyuan Zha, and Le Song. 
% "Learning Triggering Kernels for Multi-dimensional Hawkes Processes." 
% ICML (3). 2013.
%
% Luo, Dixin, et al. "Multi-task multi-dimensional hawkes processes for 
% modeling event sequences." Proceedings of the 24th International 
% Conference on Artificial Intelligence. AAAI Press, 2015.
%
% Provider:
% Hongteng Xu @ Georgia Tech
% June. 19, 2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initial 
Aest = model.A;        
muest = model.mu;
D = size(Aest, 1);

tic;
for o = 1:alg.outer
    DM = zeros(size(model.g));
    CM = DM;
    
    for n = 1:alg.inner 
        NLL = 0; % negative log-likelihood

        Amu = zeros(D, 1);
        Bmu = Amu;

        BmatA = zeros(size(Aest));
        AmatA = BmatA;
        AmatA = AmatA + 2*alg.alpha*Aest;

        % E-step: evaluate the responsibility using the current parameters    
        for c = 1:length(Seqs)
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

            Amu = Amu + Tstop - Tstart;

            dT = Tstop - Time;
            GK = Kernel_Integration_Approx(dT, model);

            Nc = length(Time);

            for i = 1:Nc

                ui = Event(i);

                AmatA(ui,:,:) = AmatA(ui,:,:)+...
                    double(Aest(ui,:,:)>0).*repmat( GK(i,:), [1,1,D] );

                ti = Time(i);             

                if n == alg.inner
                    Nums = min([ceil(dT(i)/model.dt), size(CM,1)]);
                    CM(1:Nums,:) = CM(1:Nums,:) + ...
                        repmat(sum(Aest(ui,:,:), 3), [Nums,1]);
                end
                
                lambdai = muest(ui);
                pii = muest(ui);
                pij = [];


                if i>1

                    tj = Time(1:i-1);
                    uj = Event(1:i-1);

                    dt = ti - tj;
                    gij = Kernel_Approx(dt, model);
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
                            BmatA(uuj,:,ui) = BmatA(uuj,:,ui) ...
                                + Aest(uuj,:,ui).*pij(j,:);
                            if n == alg.inner
                                Nums = min([ceil(dt(j)/model.dt), size(DM,1)]);
                                DM(Nums,:) = DM(Nums,:)+pij(j,:);
                            end
                        end
                    end
                end

                Bmu(ui) = Bmu(ui) + pii;

            end

            NLL = NLL + (Tstop - Tstart).*sum(muest);
            NLL = NLL + sum( sum( GK.*sum(Aest(Event,:,:),3) ) );

        end

        % M-step: update parameters
        mu = Bmu./Amu; 
        A = abs(sqrt(BmatA./AmatA));       
        A(isnan(A))=0;
        A(isinf(A))=0;

        Err=sum(sum(sum(abs(A-Aest))))/sum(abs(Aest(:)));
        Aest = A;
        muest = mu;
        model.A = Aest;
        model.mu = muest;
        
        fprintf('Outer=%d, Inner=%d, Obj=%f, RelErr=%f, Time=%0.2fsec\n',...
            o, n, NLL, Err, toc);

%         if Err<alg.thres || (o==alg.outer)
%             break;
%         end    
    
    end
    

        
    DM = DM./model.dt;
    CM = CM./model.g;
    gest = model.g;
    for n = 1:alg.inner_g
        for m = 1:size(model.g, 1)
            switch m
                case 1
                    a = 2*alg.alpha + CM(m,:) * (model.dt^2);
                    b = -2*alg.alpha*model.g(m+1,:);
                    c = -DM(m,:) * (model.dt^2);                    
                case size(model.g, 1)
                    a = 4*alg.alpha + CM(m,:) * (model.dt^2);
                    b = -2*alg.alpha*model.g(m-1,:);
                    c = -DM(m,:) * (model.dt^2);
                otherwise
                    a = 4*alg.alpha + CM(m,:) * (model.dt^2);
                    b = -2*alg.alpha * (model.g(m-1,:)+model.g(m+1,:));
                    c = -DM(m,:) * (model.dt^2);
            end
            
            gest(m,:) = (-b+sqrt(b.^2 - 4.*a.*c))./(2*a);
        end
        model.g = abs(gest);        
        model.g = model.g./repmat(sum(model.g),[model.M,1]);
    end                    
end


