function Loglike = Loglike_MixHP(Seqs, model, alg)


% given responsibility, calcuate the expected number of sequence belonging
% to the k-th cluster
EX = model.R;

% update parameters of Hawkes processes (mu_k, A_k), k=1,...,N
% initialize
A = model.beta;%random('exp', model.beta);
mu = sqrt(pi/2)*model.b;%random('rayl', model.b);


tmp1 = A(:)./(model.beta(:));
tmp1(isnan(tmp1)) = 0;
tmp1(isinf(tmp1)) = 0;

tmp2 = mu(:).^2./(2*model.b(:).^2);
tmp2(isnan(tmp2)) = 0;
tmp2(isinf(tmp2)) = 0;

tmp3 = log(mu(:));
tmp3(isnan(tmp3)) = 0;
tmp3(isinf(tmp3)) = 0;

Loglike = sum(tmp1)+ sum(tmp2)-sum(tmp3);





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


    for i = 1:N

        ui = Event(i);
        ti = Time(i);


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

        end

        LL = LL+log(lambdai);


    end
    LL = LL - (Tstop-Tstart).*sum(mu);
    tmp = sum(sum(repmat(G, [1,1,model.K]).*sum(A(Event,:,:,:),4),2),1);
    LL = LL - tmp(:)';


    Loglike = Loglike - EX(c,:)*LL(:);

end

Loglike = -Loglike;
