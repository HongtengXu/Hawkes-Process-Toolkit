function A = Infectivity_TVHP( T, t, Period, Shift, MaxInfect, Type )

% Generate infectivity matrix A(t) in R^{U*U}
%
% T: the time interval
% t: current time
% Period in R^{U*U}, the predefined period parameter
% Shift in R^{U*U}, the predefined time-shift parameter

switch Type
    case 1
        M = round((Period./T).*(t-Shift));
        A = MaxInfect*0.5*(1 - (-1).^M);
    case 2
        A = MaxInfect*0.5*(1+cos((2*Period./T).*t-Shift));
end

