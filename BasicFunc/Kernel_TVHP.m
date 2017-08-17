function gt = Kernel_TVHP( times, landmarks, sigma, type )

dt = repmat(landmarks(:), [1, length(times(:))])...
    - repmat(times(:)', [length(landmarks(:)), 1]);

switch type
    case 'gauss'
        gt = exp(-(dt.^2)./(2*sigma^2));
        
    case 'exp'
        gt = exp(sigma.*dt);
        gt(gt>1) = 0;
        
end