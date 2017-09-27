function Seqs = SimulationFast_Thinning_ExpHP(para, options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The fast simulation of Hawkes processes with exponential kernels
%
% Reference:
% Dassios, Angelos, and Hongbiao Zhao. 
% "Exact simulation of Hawkes process with exponentially decaying intensity." 
% Electronic Communications in Probability 18.62 (2013): 1-13.
%
% Provider:
% Hongteng Xu @ Georgia Tech
% June 13, 2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Seqs = struct('Time', [], ...
              'Mark', [], ...
              'Start', [], ...
              'Stop', [], ...
              'Feature', []);

tic
for n = 1:options.N


    t=0;
    History = [];

    lambdat = para.mu;
    mt = sum(lambdat);%SupIntensity_HP(t, History, para, options);

    while t<options.Tmax && size(History, 2)<options.Nmax

        s = random('exp', 1/mt);
        U = rand;

        lambda_ts = Intensity_Recurrent_HP(t+s, [], t, lambdat, para);
        mts = sum(lambda_ts);

        %fprintf('s=%f, v=%f\n', s, mts/mt);        
        if t+s>options.Tmax || U>mts/mt
            t = t+s;
            lambdat = lambda_ts;
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

            
            lambdat = Intensity_Recurrent_HP(t+s, index(1), t, lambdat, para);
            t = t+s;
            History=[History,[t;index(1)]];
        end
        
        mt = sum(lambdat);
        
    end
    
    Seqs(n).Time = History(1,:);
    Seqs(n).Mark = History(2,:);
    Seqs(n).Start = 0;
    Seqs(n).Stop = options.Tmax;
    
    index = find(Seqs(n).Time<=options.Tmax);
    Seqs(n).Time = Seqs(n).Time(index);
    Seqs(n).Mark = Seqs(n).Mark(index);
    
    if mod(n, 10)==0 || n==options.N
        fprintf('#seq=%d/%d, #event=%d, time=%.2fsec\n', ...
            n, options.N, length(Seqs(n).Mark), toc);
    end
end
    
    
