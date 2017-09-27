function lambda = Intensity_HP(t, History, para) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compute the intensity functions of Hawkes processes
%
% Parameters of Hawkes processes
% para.mu: base exogenous intensity
% para.A: coefficients of impact function
% para.kernel: 'exp', 'gauss'
% para.w: bandwith of kernel
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda = para.mu(:);

if ~isempty(History)
    
    Time = History(1, :);
    index = Time<=t;
    Time = Time(index);
    Event = History(2, index);
    
    
        
    basis = Kernel(t- Time(:), para);
    A = para.A(Event, :, :);
        
    for c = 1:size(para.A, 3);
        lambda(c) = lambda(c) + sum(sum(basis.*A(:,:,c)));
    end        
    
end

lambda = lambda.*double(lambda>0);
end




