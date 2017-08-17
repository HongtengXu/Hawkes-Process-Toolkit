function dis = Distance_MPP(X, Y, M, id)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Define the distance of two point processes' instances
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

dis = 0;
for n1 = 1:size(X,2)
    for n2 = 1:size(X,2)
        tmp = 1;
        for i = 1:length(id)
            tmp = tmp*( M(id(i)) - abs(X(id(i),n1)-X(id(i),n2)));
            if tmp==0
                break
            end
        end
        dis = dis + tmp;
    end
end

for n1 = 1:size(Y,2)
    for n2 = 1:size(Y,2)
        tmp = 1;
        for i = 1:length(id)
            tmp = tmp*( M(id(i)) - abs(Y(id(i),n1)-Y(id(i),n2)));
            if tmp==0
                break
            end
        end
        dis = dis + tmp;
    end
end

for n1 = 1:size(X,2)
    for n2 = 1:size(Y,2)
        tmp = 1;
        for i = 1:length(id)
            tmp = tmp*( M(id(i)) - abs(X(id(i),n1)-Y(id(i),n2)));
            if tmp==0
                break
            end
        end
        dis = dis - 2*tmp;
    end
end

dis = sqrt(abs(dis));
