function SeqsNew = SimulationFast_ConditionalThinning_ExpHP(SeqsOld, ...
                                                        para, options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The fast simulation of Hawkes processes with exponential kernels
% conditioned on history
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


SeqsNew = struct('Time', [], ...
              'Mark', [], ...
              'Start', [], ...
              'Stop', [], ...
              'Feature', []);

tic
for n = 1:length(SeqsOld)


    t=SeqsOld(n).Stop;
    History = [SeqsOld(n).Time; SeqsOld(n).Mark];

    lambdat = Intensity_HP(t, History, para);
    mt = sum(lambdat);
    
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
            History = [History,[t;index(1)]];
        end
        
        mt = sum(lambdat);
        
    end
    
    SeqsNew(n).Time = History(1,:);
    SeqsNew(n).Mark = History(2,:);
    SeqsNew(n).Start = SeqsOld(n).Stop;
    SeqsNew(n).Stop = options.Tmax;
    index = find(SeqsOld(n).Stop<=SeqsNew(n).Time & ...
        SeqsNew(n).Time<=options.Tmax);
    SeqsNew(n).Time = SeqsNew(n).Time(index);
    SeqsNew(n).Mark = SeqsNew(n).Mark(index);
    
    if mod(n, 10)==0 || n==options.N
        fprintf('#seq=%d/%d, #event=%d, time=%.2fsec\n', ...
            n, options.N, length(SeqsNew(n).Mark), toc);
    end
end
    
    
