function lambda = Intensity_TVHP(t, History, T, mu, w, Period, Shift, MaxInfect, Type) 
% INTHAWKESM  Intensity lambda(t) of a time-varying multi-dim Hawkes process
%             
% Inputs      t   - current interval 
%             History   - historical event sequences
%             T   - time interval
%             mu - U*1 vector (unconditional intensities)
%             w - parameters of decay function
%             Period, Shift - U*U parameters of infectivity matrix
% Outputs     lambda   - intensity 

if isempty(History)
    lambda = mu;
else
    Time = History(1, :);
    index = Time<=t;
    Time = Time(index);
    Event = History(2, index);

    lambda = mu;

    At = Infectivity_TVHP( T, t, Period, Shift, MaxInfect, Type );
    
    for i = 1:length(Time)
        ui = Event(i);
        
        lambda = lambda + At(:, ui).*exp(-w*(t-Time(i)));
    end   
end
