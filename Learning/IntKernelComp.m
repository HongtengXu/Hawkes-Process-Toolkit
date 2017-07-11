function KGt = IntKernelComp( tj, option, upper, lower)

dt = repmat(option.landmark(:), [1,length(tj)]) - ...
    repmat(tj(:)', [length(option.landmark),1]);



Weight = sqrt(pi/2)*option.sigmaA * ...
    exp( 0.5*option.sigmaA^2 * option.sigmaT^2 - ...
    option.sigmaT*dt);

Upp = (upper-(option.landmark(:)-option.sigmaA^2 * option.sigmaT))/(sqrt(2)*option.sigmaA);
Low = (lower-(option.landmark(:)-option.sigmaA^2 * option.sigmaT))/(sqrt(2)*option.sigmaA);

Phi = erf(Upp)-erf(Low);

KGt = Weight.*repmat(Phi, [1, length(tj)]);


    