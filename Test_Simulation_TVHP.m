%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Simulation Time-varying Hawkes Process
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

% time interval
T = 50;
% dimension
U = 2;
% intrinsic intensity
mu = [0.3;0.2];%rand(U,1)./U;
% parameter of decay function
w = 0.2;
% parameters for generating A(t)
Period = 0.3*[3,6;9,12];%2*rand(U)+1;
Shift = 0;%2*Period.^2;%.*rand(U);

Type = 2;
A = zeros(U,U,T,Type);
% the number of sequences
N = 2000;
Seqs = cell(Type,1);

for type = 1:Type%:-1:1
    if type == Type
        MaxInfect = 0.8/(U^2);
    else
        MaxInfect = 0.5/(U^2);
    end
    
    for t=1:T
        A(:,:,t,type) = Infectivity_TVHP(T, t, Period, Shift, MaxInfect, type);
    end
    figure
    for u=1:U
        for v = 1:U
            subplot(U,U,U*(u-1)+v)
            tmp=A(u,v,:,type);
            plot(1:T, tmp(:), 'r-');
            title(sprintf('a_{%d%d}(t)',u,v));
        end
    end


    Seqs{type} = Simulation_TVHP( N, T, mu, w, Period, Shift, MaxInfect, type );
    %ShowTMHP(Seq{1,type}, mu, T, w, Period, Shift, MaxInfect, type)

end
MaxInfect = 0.8/(U^2);
save('SynData_TVHP.mat','Seqs','Period','Shift','mu','w','T','MaxInfect','Type');

%%
load SynData_TVHP.mat
type = 2;
Show_Intensity_TVHP(Seqs{type}(1), mu, T, w, Period, Shift, MaxInfect, type)