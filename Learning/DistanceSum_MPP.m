function [DisS, dis] = DistanceSum_MPP(X, Y, configure)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Define the weighted distance of two point processes' instances
%
% Reference:
% Iwayama, Koji, Yoshito Hirata, and Kazuyuki Aihara. 
% "Definition of distance for nonlinear time series analysis 
% of marked point process data." 
% Physics Letters A 381.4 (2017): 257-262.
%
% Provider:
% Hongteng Xu @ Georgia Tech
% June 13, 2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dis = zeros(length(configure.weight),1);
DisS = 0;
for i = 1:length(configure.weight)
    dis(i) = Distance_MPP(X, Y, configure.M, configure.id{i});
    DisS = DisS + configure.weight(i) * dis(i);
end
    

