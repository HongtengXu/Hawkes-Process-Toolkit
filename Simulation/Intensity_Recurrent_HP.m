function lambda = Intensity_Recurrent_HP(t_current, event_current, t_old, lambdat, para) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compute the intensity functions of Hawkes processes recurrently according
% to historical intensity and current event.
%
% Parameters of Hawkes processes
% para.mu: base exogenous intensity
% para.A: coefficients of impact function
% para.kernel: 'exp', 'gauss'
% para.w: bandwith of kernel
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt = t_current - t_old;
weight = exp(-para.w * dt);
lambda = para.mu + (lambdat - para.mu).*weight;


if ~isempty(event_current)
    
    lambda = lambda + reshape(para.A(event_current,1,:), [size(para.A,3),1]);
        
end

end




