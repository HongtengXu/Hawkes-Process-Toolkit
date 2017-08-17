function Loglike = Loglike_HP_ODE( Seqs, model, alg )
                                                        
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

tic;


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



        ti = Time(i);             



        lambdai = muest(ui);



        if i>1

            tj = Time(1:i-1);
            uj = Event(1:i-1);

            dt = ti - tj;
            gij = Kernel_Approx(dt, model);
            auiuj = Aest(uj, :, ui);
            pij = auiuj .* gij;
            lambdai = lambdai + sum(pij(:));
        end

        Loglike = Loglike - log(lambdai);




    end

    Loglike = Loglike + (Tstop - Tstart).*sum(muest);
    Loglike = Loglike + sum( sum( GK.*sum(Aest(Event,:,:),3) ) );

end

Loglike = -Loglike;
        
    

