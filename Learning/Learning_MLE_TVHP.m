function model = Learning_MLE_TVHP( Seqs, option )

tic
% initialize
model.A = rand(option.D, option.M, option.D);
model.mu = rand(option.D, 1);

Z = model.A;
U = zeros(size(Z));

for o = 1:option.outer
    rho = option.rho * option.factor^o;
    
    % update model.A, mu
    for l = 1:option.inner
        
        NLLinner = 0;
        AA = rho*ones(size(model.A));
        AB = rho*(U-Z);
        AC = zeros(size(model.A));
                
        muD = 0;
        muN = zeros(size(model.mu));
        
        Amc = sum(model.A,3);
        for n = 1:length(Seqs)
            
            Times = Seqs(n).Time;
            Tstart = Seqs(n).Start;
            Tstop = Seqs(n).Stop;
            Events = Seqs(n).Mark;

            
            muD = muD + (Tstop-Tstart);
            
            In = length(Times);
            kappat = Kernel_TVHP(Times, option.landmark, ...
                option.sigmaA, option.typeA);
            
            
            for i = 1:In
                ti = Times(i);
                ci = Events(i);
                lambdai = model.mu(ci);
                pii = model.mu(ci);
                
                
                lower = ti;
                if i==In
                    upper = Tstop;
                else
                    upper = Times(i+1);
                end
                KGt = IntKernelComp_TVHP( Times(1:i), option, upper, lower);
                for j = 1:i                        
                    AB(Events(j), :, :) = AB(Events(j), :, :) +...
                        repmat(KGt(:,j)', [1,1,option.D]);
                end
                NLLinner = NLLinner + sum(sum(KGt'.*Amc(Events(1:i),:)));
                
                if i>1
                    tj = Times(1:i-1);
                    cj = Events(1:i-1);
                    
                    acicj = model.A(cj,:,ci);
                    gt = Kernel_TVHP(ti, tj, option.sigmaT, option.typeT);                
                    
                    pijm = repmat(gt(:), [1,option.M])...
                        .* acicj...
                        .* repmat(kappat(:,i)', [i-1, 1]);
                    
                    lambdai = lambdai + sum(pijm(:));
                
                    pijm = pijm./lambdai;
                    
                    for j = 1:i-1
                        AC(cj(j), :, ci) = AC(cj(j), :, ci) - pijm(j, :);                        
                    end
                end
                
                NLLinner = NLLinner - log(lambdai);
                
                pii = pii./lambdai;
                muN(ci) = muN(ci)+pii;
                
                
            end
            
            NLLinner = NLLinner + (Tstop-Tstart)*sum(model.mu(:));
            
        end
        
        
        tmpA = (-AB+sqrt(AB.^2 - 4*AA.*AC))./(2*AA);       
        tmpmu = muN./muD;
        
        RelErr = norm(tmpA(:)-model.A(:))/(option.D^2*option.M);
        
        fprintf('Update A, mu: Outer%d, Inner%d, NLL=%f, NLL+Reg=%f, Err=%d, Time=%.2fsec\n', ...
            o, l, NLLinner, NLLinner+option.alphaS*sum(abs(model.A(:))), RelErr, toc);
        
        
        if RelErr<option.thres
            model.A = tmpA;
            model.mu = tmpmu;
            break
        else
            model.A = tmpA;
            model.mu = tmpmu;
        end
        
        
    end
    
    % update Z
    Z = SoftThreshold_S( model.A+U, option.alphaS/rho );
    
    % update U
    U = U + (model.A - Z);
    
end