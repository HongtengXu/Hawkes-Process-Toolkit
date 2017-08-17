function Z=SoftThreshold_GS( A, thres )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Soft-thresholding for group lasso
%
% Reference:
% Yuan, Ming, and Yi Lin. 
% "Model selection and estimation in regression with grouped variables." 
% Journal of the Royal Statistical Society: Series B 
% (Statistical Methodology) 68.1 (2006): 49-67.
%
% Provider:
% Hongteng Xu @ Georgia Tech
% June 12, 2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Z = zeros(size(A));
for u = 1:size(A, 3)
    for v = 1:size(A, 1)
        tmp = 1 - thres/norm(A(v,:,u));
        if tmp>0
            Z(v,:,u) = tmp*A(v,:,u);
        end
    end
end
