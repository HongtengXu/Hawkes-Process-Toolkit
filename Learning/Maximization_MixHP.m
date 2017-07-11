function [model, NLL] = Maximization_MixHP(Seqs, model, alg)


% given responsibility, calcuate the expected number of sequence belonging
% to the k-th cluster
EX = model.R;

% update parameters of Hawkes processes (mu_k, A_k), k=1,...,N
% initialize
A = model.beta;%random('exp', model.beta);
mu = sqrt(pi/2)*model.b;%random('rayl', model.b);


for in = 1:alg.inner
    tmp1 = A(:)./(model.beta(:));
    tmp1(isnan(tmp1)) = 0;
    tmp1(isinf(tmp1)) = 0;

    tmp2 = mu(:).^2./(2*model.b(:).^2);
    tmp2(isnan(tmp2)) = 0;
    tmp2(isinf(tmp2)) = 0;

    tmp3 = log(mu(:));
    tmp3(isnan(tmp3)) = 0;
    tmp3(isinf(tmp3)) = 0;

    NLL = sum(tmp1)+ sum(tmp2)-sum(tmp3);


    MuA = 1./(model.b.^2);
    MuA(isinf(MuA)) = 0;
    MuB = 0;
    MuC = -ones(size(model.b));


    AB = zeros(size(A));
    AA = 1./(model.beta);
    AA(isinf(AA)) = 0;


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

        N = length(Time);
        % calculate the integral decay function in the log-likelihood function
        G = Kernel_Integration(Tend - Time, model);


        TMPAA = zeros(size(A));
        TMPAB = zeros(size(A));
        %TMPMuB = zeros(size(mu));
        TMPMuC = zeros(size(mu));
        LL = 0;
            
        for i = 1:N

            ui = Event(i);
            ti = Time(i);
            TMPAA(ui,:,:,:) = TMPAA(ui,:,:,:)+ ...
                        repmat(G(i,:), [1, 1, model.K, model.D]);


            lambdai = mu(ui,:)+eps;
            pii = lambdai;
                
            if i>1
                tj = Time(1:i-1);
                uj = Event(1:i-1);

                gij = Kernel(ti-tj, model);
                auiuj = A(uj, :, :, ui);
                pij = repmat(gij, [1,1,model.K,1]).* auiuj;
                    
                    
                tmp = sum(sum(pij,1),2);
                lambdai = lambdai + tmp(:)';
                    
                    
                pij = pij./repmat(reshape(lambdai,[1,1,model.K]),...
                        [size(pij,1),size(pij,2), 1]);
                    
                for j=1:i-1
                    uj = Event(j);
                    TMPAB(uj,:,:,ui) = TMPAB(uj,:,:,ui) - pij(j,:,:);
                end
                 
            end

            LL = LL+log(lambdai);
                
                
                
            pii = pii./lambdai;                
            TMPMuC(ui,:)=TMPMuC(ui,:)-pii;
                
        end
        LL = LL - Tend.*sum(mu);
        tmp = sum(sum(repmat(G, [1,1,model.K]).*sum(A(Event,:,:,:),4),2),1);
        LL = LL - tmp(:)';

        % XX = (LL - max(LL));
        % EX(c,:)=(model.p'.*(exp(XX)+options.bias))./((exp(XX)+options.bias)*model.p);
            
        MuB = MuB + Tend*EX(c,:);
        for k=1:model.K
            AA(:,:,k,:) = AA(:,:,k,:) + EX(c,k)*TMPAA(:,:,k,:);
            AB(:,:,k,:) = AB(:,:,k,:) + EX(c,k)*TMPAB(:,:,k,:);
            MuC(:,k) = MuC(:,k) + EX(c,k)*TMPMuC(:,k);                
        end
        NLL = NLL - EX(c,:)*LL(:);
            
    end
        
    MuBB = repmat(MuB,[model.D, 1]);    
    % M-step: update parameters
    mutmp = (-MuBB+sqrt(MuBB.^2 - 4*MuA.*MuC))./(2*MuA);
    Atmp = -AB./AA;  
        
        
    Atmp(isnan(Atmp)) = 0;
    Atmp(isinf(Atmp)) = 0;
    mutmp(isnan(mutmp)) = 0;
    mutmp(isinf(mutmp)) = 0;
        
    % check convergence
    Err=sum(abs(A(:)-Atmp(:)))/sum(abs(A(:)));
    fprintf('Inner=%d, Obj=%f, RelErr=%f, Time=%0.2fsec\n',...
            in, NLL, Err, toc);

    A = Atmp;
    mu = mutmp;
    if Err<alg.thres || in==alg.inner
        break;
    end    
    
end    

model.beta = A;
model.b = sqrt(2/pi)*mu;

