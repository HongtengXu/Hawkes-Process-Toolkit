function mt = SupIntensity_HP(t, History, para, options) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compute the super bound of intensity function of Hawkes processes
%
% Parameters of Hawkes processes
% para.mu: base exogenous intensity
% para.A: coefficients of impact function
% para.kernel: 'exp', 'gauss'
% para.w: bandwith of kernel
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if isempty(History)
    mt = sum(para.mu);
else
    Time = History(1, :);
    index = Time<=t;
    Time = Time(index);
    Event = History(2, index);
    
    MT = sum(para.mu)*ones(1, options.M);
    for m=1:options.M
        t_current = t+(m-1)*options.tstep/options.M;
        
        basis = Kernel(t_current-Time(:), para);
        A = para.A(Event, :, :);
        
        for c = 1:size(para.A, 3);
            MT(m) = MT(m) + sum(sum(basis.*A(:,:,c)));
        end        
    end  
    mt = max(MT);
end

end




