function [A, Phi] = ImpactFunc_ODE( para )


Phi = zeros(size(para.A, 1), para.M, size(para.A, 3));
A = zeros(size(para.A, 3), size(para.A, 1));


for u = 1:size(para.A, 3)
    for v = 1:size(para.A, 1)
        
        Phi(v,:,u) = para.A(v,:,u)*para.g';
        A(u, v) = sum(Phi(v,:,u))*para.dt;
        
    end
end