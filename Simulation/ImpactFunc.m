function [A, Phi] = ImpactFunc( para, options )


Phi = zeros(size(para.A, 1), options.M, size(para.A, 3));
A = zeros(size(para.A, 3), size(para.A, 1));

time_stamp = 0:options.dt:(options.M-1)*options.dt;

for u = 1:size(para.A, 3)
    for v = 1:size(para.A, 1)
        basis_int = Kernel_Integration(options.Tmax, para);
        A(u, v) = para.A(v,:,u)*basis_int(:);
        
        basis = Kernel(time_stamp(:), para);
        Phi(v,:,u) = para.A(v,:,u)*basis';
        
    end
end