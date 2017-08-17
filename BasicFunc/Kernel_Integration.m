function G = Kernel_Integration(dt, para)


% dt = t_current - t_hist(:);
distance = repmat(dt(:), [1, length(para.landmark(:))]) - ...
            repmat(para.landmark(:)', [length(dt), 1]);
landmark = repmat(para.landmark(:)', [length(dt), 1]);

switch para.kernel
    case 'exp'
        G = 1 - exp(-para.w * (distance-landmark));
        G(G<0) = 0;
        
    case 'gauss'
        G = 0.5*( erf(distance./(sqrt(2)*para.w)) ...
            + erf(landmark./(sqrt(2)*para.w)) );
        
    otherwise
        disp('Error: please assign a kernel function!');
end
    
