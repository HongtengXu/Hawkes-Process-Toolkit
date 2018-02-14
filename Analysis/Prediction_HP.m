function [lambda_current, count_expect] = Prediction_HP(t_start, t_end, ...
                                                    History, para, options)
                                                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Prediction of current intensity function and expected number of future 
% events.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda_current = Intensity_HP(t_start, History, para);
count_expect = zeros(size(para.A, 1), 1);


for n = 1:options.NumTest
    Hist = History;
    t = t_start;
    
    mt = SupIntensity_HP(t, Hist, para, options);

    while t<t_end && size(Hist, 2)<options.Nmax

        s = random('exp', 1/mt);
        U = rand;

        lambda_ts = Intensity_HP(t+s, Hist, para);
        mts = sum(lambda_ts);

        %fprintf('s=%f, v=%f\n', s, mts/mt);        
        if t+s>t_end || U>mts/mt
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
            Hist=[Hist,[t;index(1)]];
            count_expect(index(1)) = count_expect(index(1)) + 1;
        end
        
        mt = SupIntensity_HP(t, Hist, para, options);
        
    end
    
    if mod(n, 10)==0 || n==options.NumTest
        fprintf('#seq=%d/%d, #event=%d, time=%.2fsec\n', ...
            n, options.NumTest, size(Hist,2), toc);
    end
end

count_expect = count_expect./options.NumTest;

