function Loglike = Loglike_Basis( Seqs, model, alg )
                                                        
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



D = size(Aest, 1);


tic;

        
Loglike = 0; % negative log-likelihood



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
    GK = Kernel_Integration(dT, model);

    Nc = length(Time);

    for i = 1:Nc

        ui = Event(i);


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

        Loglike = Loglike - log(lambdai);

    end

    Loglike = Loglike + (Tstop-Tstart).*sum(muest);
    Loglike = Loglike + sum( sum( GK.*sum(Aest(Event,:,:),3) ) );



end

Loglike = -Loglike;
                
        