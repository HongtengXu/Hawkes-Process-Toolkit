function model = Learning_LS_Discrete( Seqs, model, alg )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Learning Hawkes processes via least squares-based nonparametric method.
%
% Reference:
% Eichler, Michael, Rainer Dahlhaus, and Johannes Dueck. 
% "Graphical modeling for multivariate hawkes processes 
% with nonparametric link functions." 
% Journal of Time Series Analysis (2016).
%
% Provider:
% Hongteng Xu @ Georgia Tech
% June. 12, 2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


N = length(Seqs);
dt = model.h/4;
Mu = 0;
Phi = 0;

for n = 1:4
    Gamma_h = zeros(model.D * model.k);
    gamma_h = zeros(model.D, model.D*model.k);
    YH = zeros(model.D, 1);
    YHK = zeros(model.D*model.k, 1);

    for i = 1:length(Seqs)
        Time = Seqs(i).Time;
        Event = Seqs(i).Mark;
        Tstart = Seqs(i).Start;
            
        if isempty(model.Tmax)
            Tstop = Seqs(i).Stop;
        else
            Tstop = model.Tmax;
            indt = Time < model.Tmax;
            Time = Time(indt);
            Event = Event(indt);
        end
            
        Th = floor((Tstop-Tstart) ./ model.h);
        Y = zeros(model.D, Th);
        for d = 1:model.D
            ind = Event == d;
            Time_d = Time(ind) - Tstart;
            for t = 1:Th
                td = find(Time_d>=(t-1)*model.h + (n-1)*dt & Time_d<t*model.h + (n-1)*dt);
                Y(d, t) = length(td);
            end
        end
        Thk = Th - model.k;

        if Thk>0
            Yk = zeros(model.D * model.k, Th-model.k);
            for k = 1:Th-model.k
                Yk(:, k) = reshape(Y(:,k+model.k-1:-1:k), [model.D * model.k, 1]);
            end


            Yh = mean(Y(:, model.k+1:Th), 2);
            Yhk = mean(Yk, 2);
            YH = YH + (1/N)*Yh;
            YHK = YHK + (1/N)*Yhk;



            Gamma_tmp = (1/(N*Thk))*(Yk - repmat(Yhk,[1,Thk]))*...
                        (Yk - repmat(Yhk, [1,Thk]))';
            Gamma_h = Gamma_h + Gamma_tmp; 

            gamma_tmp = (1/(N*Thk))*(Y(:, model.k+1:Th) - repmat(Yh,[1,Thk]))*...
                        (Yk - repmat(Yhk, [1,Thk]))';
            gamma_h = gamma_h + gamma_tmp;
        end
    end


    phi = (gamma_h / Gamma_h)./model.h;
    mu = YH - phi * YHK;

    Mu = Mu + mu;
    Phi = Phi + phi;
end

Mu = Mu/4;
Phi = Phi/4;

model.mu = Mu./model.h;
model.A = zeros(model.D, model.k, model.D);
for d = 1:model.D
    for k = 1:model.k
        model.A(:, k, d) = Phi(d, 1+(k-1)*model.D:k*model.D)';
    end
end

model.GCG = reshape(sum(model.A, 2), [model.D, model.D])';
% end
% 
% function f = objfun(x, A, B)
% 
% x = reshape(x, [size(B,1), size(A,1)]);
% f = norm(x*A-B, 'fro')^2;
% 
% end