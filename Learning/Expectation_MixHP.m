function model = Expectation_MixHP(Seqs, model, alg)

Nk = sum(model.R,1); % 10.51
alpha = model.alpha+Nk; % Dirichlet
Elogpi = psi(0,alpha)-psi(0,sum(alpha)); % 10.66

    
% calculate responsibility
tic;   
EX = zeros(length(Seqs), model.K);
       

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
            
    N = length(Time);
    % calculate the integral decay function in the log-likelihood function
    G = Kernel_Integration(Tstop - Time, model);


    LL = Elogpi;

    for i = 1:N

        ui = Event(i);
        ti = Time(i);
        
        Elambdai = sqrt(pi/2).*model.b(ui,:)+eps;
        Vlambdai = (2-pi/2).*(model.b(ui,:)).^2;
        if i>1
            tj = Time(1:i-1);
            uj = Event(1:i-1);

            gij = Kernel(ti-tj, model);
            auiuj = model.beta(uj, :, :, ui);
            pij = repmat(gij, [1,1,model.K,1]).* auiuj;

            tmp = sum(sum(pij,1),2);
            Elambdai = Elambdai + tmp(:)';  
            tmp = sum(sum(pij.^2, 1), 2);
            Vlambdai = Vlambdai + tmp(:)';
        end

        LL = LL+log(Elambdai) - Vlambdai./(2*Elambdai.^2);




    end
    LL = LL - (Tstop-Tstart).*sqrt(pi/2).*sum(model.b);
    tmp = sum(sum(repmat(G, [1,1,model.K]).*...
        sum(model.beta(Event,:,:,:),4),1),2);
    LL = LL- tmp(:)';

    XX = (LL - max(LL));
    EX(c,:)=(exp(XX))./sum(exp(XX));

end

model.R = EX;

end