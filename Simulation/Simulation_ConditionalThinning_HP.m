function SeqsNew = Simulation_ConditionalThinning_HP(...
                                                SeqsOld, para, options)
                                            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Implementation of Ogata's thinning method to simulate Hawkes processes
% conditioned on Historical event sequences (SeqsOld)
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


SeqsNew = struct('Time', [], ...
              'Mark', [], ...
              'Start', [], ...
              'Stop', [], ...
              'Feature', []);
tic
for n = 1:length(SeqsOld)


    t=SeqsOld(n).Stop;
    History = [SeqsOld(n).Time; SeqsOld(n).Mark];

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
    %History
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