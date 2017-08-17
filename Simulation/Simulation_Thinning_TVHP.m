function Seqs = Simulation_Thinning_TVHP( para, options )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Simulate event sequence of time-varying multi-dimensional Hawkes process
%
% N: the number of sequences
% T: the time interval of each sequence
% mu: intrinsic intensity vector
% w: parameter of decay function
% Period, Shift: parameters of infectivity matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Seqs = struct('Time', [], ...
              'Mark', [], ...
              'Start', [], ...
              'Stop', [], ...
              'Feature', []);

tic
for ns = 1:options.N


    t=0;
    History = [];
    
    lambda_t = Intensity_TVHP( t, History, para, options );
    mt = sum(lambda_t);
    
    while t<options.Tmax && size(History, 2)<options.Nmax
        
        s = random('exp', 1/mt);
        U = rand;

        
        lambda_ts = Intensity_TVHP( t+s, History, para, options );
        mts = sum(lambda_ts);


        %fprintf('s=%f, v=%f\n', s, mts/mt);        
        if t+s>options.Tmax || U>mts/mt
            t = t+s;
        else
            u = rand*mts;
            sumIs = 0;
            for d=1:length(lambda_ts)
                sumIs = sumIs + lambda_ts(d);
                if sumIs >= u
                    break;
                end
            end
            index = d;

            t = t+s;
            History=[History,[t;index(1)]];
        end
        
        lambda_t = Intensity_TVHP(t, History, para, options);
        mt = sum(lambda_t);
    end
    
    Seqs(ns).Time = History(1,:);
    Seqs(ns).Mark = History(2,:);
    Seqs(ns).Start = 0;
    Seqs(ns).Stop = options.Tmax;
    
    if mod(ns, 10)==0 || ns==options.N
        fprintf('#seq=%d/%d, #event=%d, time=%.2fsec\n', ...
            ns, options.N, size(History,2), toc);
    end
end