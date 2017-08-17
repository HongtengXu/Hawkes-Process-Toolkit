function Seqs = Simulation_TVHP( N, T, mu, w, Period, Shift, MaxInfect, Type )

% Simulate event sequence of time-varying multi-dimensional Hawkes process
%
% N: the number of sequences
% T: the time interval of each sequence
% mu: intrinsic intensity vector
% w: parameter of decay function
% Period, Shift: parameters of infectivity matrix

Seqs = struct('Time', [], ...
              'Mark', [], ...
              'Start', [], ...
              'Stop', [], ...
              'Feature', []);

tic
for ns=1:N
%     U = length(mu);
%     h = -log(rand(U,1))./mu;
%     t = max(h);

    t=0;
    History = [];

    lambda_t = Intensity_TVHP(t, History, T, mu, w, Period, Shift, MaxInfect, Type);
    mt = sum(lambda_t);
    
    while t<T
        
        s = random('exp', 1/mt);
        U = rand;

        lambda_ts = Intensity_TVHP(t+s, History, T, mu, w, Period, Shift, MaxInfect, Type);
        mts = sum(lambda_ts);


        if t+s>T || U>mts/mt
            
            t = t+s;
            mt = mts;
            
        else

            u = rand()*mts;
            sumIs = 0;
            for d=1:length(lambda_ts)
               sumIs = sumIs + lambda_ts(d);
               if sumIs >= u
                   break;
               end
            end
            index = d;
%             prob=lambda_ts./mts;
%             pd = makedist('Multinomial','probabilities', prob');
%             index = random(pd,1);

            t = t+s;
            History=[History,[t;index(1)]];

            lambda_t = Intensity_TVHP(t, History, T, mu, w, Period, Shift, MaxInfect, Type);
            mt = sum(lambda_t);

        end
    end
    Seqs(ns).Time = History(1,:);
    Seqs(ns).Mark = History(2,:);
    Seqs(ns).Start = 0;
    Seqs(ns).Stop = T;
    
    if mod(ns, 10)==0 || ns==N
        fprintf('Type=%d, #seq=%d/%d, #event=%d, time=%.2fsec\n', ...
            Type, ns, N, size(History,2), toc);
    end
end