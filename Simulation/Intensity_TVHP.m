function lambda = Intensity_TVHP( time, history, model, options )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% t: current time
% mu: base intensity in R^U
% At: time varying infectivity A(t) in R^{U*U}
% history: historical events
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



lambda = repmat( model.mu(:), [1, length(time)] );
if ~isempty(history)
    
    switch options.type
        case 'Est'
            At = Infectivity_TV_Est(time, model);
        case 'Syn'
            At = Infectivity_TV_Syn(T, time, model);      
    end
    
    for n = 1:length(time)

        index = find(history(1,:)<=time(n));
        decay = options.w*exp(-options.w*(time(n)-history(1,index)));
        eventID = history(2, index);
        
        for u = 1:size(At,1)
            lambda(u, n) = lambda(u, n) + decay(:)'*At(eventID, n, u);
        end
    end

end

