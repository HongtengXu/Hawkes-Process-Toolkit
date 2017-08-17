function phi = ImpactFunction( u, dt, para )

A = reshape(para.A(u,:,:), [size(para.A, 2), size(para.A, 3)]);
basis = Kernel( dt, para );
phi = A'*basis(:);

