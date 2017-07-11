function model = Learning_Cluster_Basis(Seqs, model, alg)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The Dirichlet Process Hawkes Mixture Processes via Variational Inference.
% 
% The unified framework is:
% The intensity of the u-th event type in the k-th process at time t:
%  lambda_{u,k}(t)=mu_{u,k}+sum_{t_j<t}sum_{m=1}^{M}a_{u,u_j,m,k}G_m(t-t_j)
%
% The joint distribution of observed sequence and unknown parameters is
% p(X, Z, w, mu, A) = p(X | Z, mu, A) p(Z | w) p(w) p(mu) p(A)
%
% p(X | Z, mu, A) = prod_n prod_k HP(x_n | mu_k, A_k)^{z_{nk}}
% p(Z | w) = prod_n prod_k w_k^{z_{nk}}
% p(w) = Dir(w | alpha_0) = C(alpha_0) prod_k w_k^{alpha0-1}
% p(mu) = prod_k Rayleigh( mu_k | b_k )
% p(A) = prod_k prod_u prod_u prod_d exp(a_{kuu'd} | v)
%
% We learn the model via variantional inference
% Given variantional distribution: 
%   q(Z, w, mu, A) = q(Z)q(w, mu, A) = q(Z)q(w)prod_k q(mu_k)q(A_k)
% E-step:
%   max ln q(Z) => max E_{w, mu, A}[ln p(X, Z, w, mu, A)]
%               => max E_{w}[ln p(Z | w)] + E_{mu, A}[ln p(X|Z, mu, A)]
%
% M-step: 
%   max ln q(w, mu, A) 
% => max ln p(w) + sum_k ln p(mu, A) + E_{Z}[ln p(Z|w)]
%        + sum_k sum_n E[z_{nk}] ln HP(x_n | mu_k, A_k)
% 
% 
% Reference:
% Xu, Hongteng, and Hongyuan Zha. 
% "A Dirichlet Mixture Model of Hawkes Processes 
% for Event Sequence Clustering." 
% arXiv preprint arXiv:1701.09177 (2017). 
%
% Provider:
% Hongteng Xu @ Georgia Tech
% June 12, 2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;

NLL = zeros(alg.outer, 1);
for o = 1:alg.outer
    [model, NLL(o)] = Maximization_MixHP(Seqs, model, alg);
    model = Expectation_MixHP(Seqs, model, alg);
    
    fprintf('MixMHP: Iter=%d, Obj=%f, Time=%0.2fsec\n', o, NLL(o), toc)  
end
model.NLL = NLL;

