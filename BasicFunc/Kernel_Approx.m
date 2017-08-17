function g = Kernel_Approx(dt, para)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compute the value of kernel function at different time
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


g = zeros(length(dt(:)), size(para.g,2));

M = size(para.g,1);
Nums = ceil(dt./para.dt);
for i = 1:length(dt(:))
    if Nums(i)<=M
        g(i,:) = para.g(Nums(i),:);
    else
        g(i,:) = 0;
    end
end


    
