function At = Infectivity_TV_Syn(T, t, para)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Generate infectivity matrix A(t) in R^{U*U}
%
% T: the time interval
% t: current time
% para.Period in R^{1*U*U}, the predefined period parameter
% para.Shift in R^{1*U*U}, the predefined time-shift parameter
% para.max_infect, the maximum infectivity
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L = length(t);
period = repmat(para.period, [1, L, 1]);
shift = repmat(para.shift, [1, L, 1]);
tt = repmat(t(:)', [size(period,1), 1, size(period,3)]);

switch para.type
    case 'gate'
        M = round((period./T).*(tt - shift));
        At = para.max_infect*0.5*(1 - (-1).^M);
    case 'sine'
        At = para.max_infect*0.5*(1+sin((2*period./T).*(t - shift)));
end

