function History = Simulation_Thinning_Poisson(mu, t_start, t_end)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Implement thinning method to simulate homogeneous Poisson processes
%
% Provider:
% Hongteng Xu @ Georgia Tech
% June 10, 2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tic

t=t_start;
History = [];

mt = sum(mu);

while t<t_end
    s = random('exp', 1/mt);
    t = t+s;

    u = rand*mt;
    sumIs = 0;
    for d=1:length(mu)
        sumIs = sumIs + mu(d);
        if sumIs >= u
            break;
        end
    end
    index = d;

    History=[History,[t;index(1)]];

end


    
    
