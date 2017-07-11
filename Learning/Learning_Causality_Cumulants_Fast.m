function model = Learning_Causality_Cumulants_Fast(Seqs, model, alg)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Learning Granger causality graph based on cumulants of Hawkes processes
%
% Reference:
% Achab, Massil, et al. 
% "Uncovering Causality from Multivariate Hawkes Integrated Cumulants." 
% arXiv preprint arXiv:1607.06333 (2016).
%
% Provider:
% Hongteng Xu @ Georgia Tech
% June 16, 2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tic;
disp('Calculate cumulants...')
model = Estimate_Cumulants_2(Seqs, model);
fprintf('Time=%.2fsec\n', toc);

%R = randn(model.D);
L = diag(model.Cumulant1);
I = eye(model.D);
Lhat = model.Cumulant1;
Chat = model.Cumulant2;

L2 = diag(model.Cumulant1.^(-0.5));
[U,S,V] = svd(Chat);
C2 = U*(S.^0.5)*V';
R = C2 * L2;

V = repmat(-2*Lhat', [model.D, 1]);

%alpha1 = 1/norm(Chat, 'fro')^2;

C = R*L*R';
K = (R.*R)*Chat' + 2*(R.*(Chat-R*L))*R';
fR = alpha2 * norm(K - Khat, 'fro')^2 + ...
     alpha1 * norm(C - Chat, 'fro')^2;

Lreal = diag(model.Rreal*model.mureal);
Creal = model.Rreal*Lreal*model.Rreal';
Kreal = (model.Rreal.*model.Rreal)*Creal' + ...
    2*(model.Rreal.*(Creal-model.Rreal*Lreal))*model.Rreal'; 
freal = alpha2 * norm(Kreal - Khat, 'fro')^2 + ...
     alpha1 * norm(Creal - Chat, 'fro')^2;
fprintf('Iter=%d, Obj=%.4f, Obj_min=%.4f, Time=%.2fsec\n', 0, fR, freal, toc);
for nn = 1:alg.iterNum
    lr = alg.lr * (0.999)^(nn-1);
    
    pCpR = alpha1*4*L*(R*R')*L*R;
    pKpR = zeros(model.D);
    
    for m = 1:model.D
        for n = 1:model.D
            tmp = zeros(model.D);
            for i = 1:model.D
                for j = 1:model.D
                    if m==i && m~=j
                        tmp(i,j) = 2*R(m,n)*(Chat(j,n)+R(j,n)*V(j,n)) + ...
                            2*Chat(m,n)*R(j,n);
                    end
                    
                    if m==j && m~=i
                        tmp(i,j) = R(i,n)^2 * V(i,n) + 2*Chat(i,n)*R(i,n);
                    end
                    
                    if m==i && m==j
                        tmp(i,j) = 2*R(m,n)*C(m,n) + 3*R(m,n)^2 * V(m,n) + ...
                            4*Chat(m,n)*R(m,n);
                    end
                end 
            end
            
            pKpR(m,n) = trace(2*(K-Khat)' * tmp);
            
        end
    end
    pKpR = alpha2*pKpR;
    
    Rnew = R - lr * (pKpR + pCpR);
    
    err = norm(Rnew - R, 'fro')/norm(R, 'fro');
    R = Rnew;
    
    C = R*L*R';
    K = (R.*R)*(Chat + R.*V)' + (C.*R)*R';
    fR = alpha2 * norm(K - Khat, 'fro')^2 + ...
            alpha1 * norm(C - Chat, 'fro')^2;
        
    fprintf('Iter=%d, Obj=%.4f, Obj_min=%.4f, Err=%.4f, Time=%.2fsec\n', ...
        nn, fR, freal, err, toc);
end

disp('Get Granger causality graph via I-inverse(R)...')
model.R = R;
model.GCG = I - inv(R');
fprintf('Finish! Time=%.2fsec\n', toc);


