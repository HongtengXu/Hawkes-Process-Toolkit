function Seqs = Simulation_Branch_HP(para, options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Simulate Hawkes processes as Branch processes
%
% Reference:
% Møller, Jesper, and Jakob G. Rasmussen. 
% "Approximate simulation of Hawkes processes." 
% Methodology and Computing in Applied Probability 8.1 (2006): 53-64.
%
% Provider:
% Hongteng Xu @ Georgia Tech
% June 10, 2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Seqs = struct('Time', [], ...
              'Mark', [], ...
              'Start', [], ...
              'Stop', [], ...
              'Feature', []);
tic
for n = 1:options.N
    
    % the 0-th generation, simulate exogeneous events via Poisson processes
    History = Simulation_Thinning_Poisson(para.mu, 0, options.Tmax);
    current_set = History;
    
    for k = 1:options.GenerationNum
        future_set = [];
        for i = 1:size(current_set, 2)
            ti = current_set(1,i);
            ui = current_set(2,i);
            t = 0;
            
            phi_t = ImpactFunction(ui, t, para);
            mt = sum(phi_t);
            
            while t<options.Tmax-ti

                s = random('exp', 1/mt);
                U = rand;

                phi_ts = ImpactFunction(ui, t+s, para);
                mts = sum(phi_ts);

                %fprintf('s=%f, v=%f\n', s, mts/mt);        
                if t+s>options.Tmax-ti || U>mts/mt
                    t = t+s;
                else
                    u = rand*mts;
                    sumIs = 0;
                    for d=1:length(phi_ts)
                        sumIs = sumIs + phi_ts(d);
                        if sumIs >= u
                            break;
                        end
                    end
                    index = d;

                    t = t+s;
                    future_set=[future_set,[t+ti;index(1)]];
                end
        
                phi_t = ImpactFunction(ui, t, para);
                mt = sum(phi_t);
            end            
        end
        
        if isempty(future_set)
            break
        else
            current_set = future_set;
            History = [History, current_set];
        end
    end
    
    [~, index] = sort(History(1,:), 'ascend');
    Seqs(n).Time = History(1,index);
    Seqs(n).Mark = History(2,index);
    Seqs(n).Start = 0;
    Seqs(n).Stop = options.Tmax;
    
    if mod(n, 10)==0 || n==options.N
        fprintf('#seq=%d/%d, #event=%d, time=%.2fsec\n', ...
            n, options.N, size(History,2), toc);
    end
end
    
    
