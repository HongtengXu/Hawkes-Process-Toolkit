function Seqs = Simulation_Thinning_HP(para, options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Implementation of Ogata's thinning method to simulate Hawkes processes
%
% Reference:
% Ogata, Yosihiko. "On Lewis' simulation method for point processes." 
% IEEE Transactions on Information Theory 27.1 (1981): 23-31.
%
% Provider:
% Hongteng Xu @ Georgia Tech
% June 10, 2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Seqs = struct('Time', [], 'Mark', []);

tic
for n = 1:options.N


    t=0;
    History = [];

    mt = SupIntensity_HP(t, History, para, options);

    while t<options.Tmax && size(History, 2)<options.Nmax

        s = random('exp', 1/mt);
        U = rand;

        lambda_ts = Intensity_HP(t+s, History, para);
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
        
        mt = SupIntensity_HP(t, History, para, options);
        
    end
    Seqs(n).Time = History(1,:);
    Seqs(n).Mark = History(2,:);

    if mod(n, 10)==0 || n==options.N
        fprintf('#seq=%d/%d, #event=%d, time=%.2fsec\n', ...
            n, options.N, size(History,2), toc);
    end
end
    
    
