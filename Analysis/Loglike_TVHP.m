function Loglike = Loglike_TVHP( Seqs, model, option )


Loglike = 0;     
Amc = sum(model.A,3);

for n = 1:length(Seqs)

    Times = Seqs(n).Time;
    Tstart = Seqs(n).Start;
    Tstop = Seqs(n).Stop;
    Events = Seqs(n).Mark;

    In = length(Times);
    kappat = Kernel_TVHP(Times, option.landmark, option.sigmaA, option.typeA);


    for i = 1:In
        ti = Times(i);
        ci = Events(i);
        lambdai = model.mu(ci);


        lower = ti;
        if i==In
            upper = Tstop;
        else
            upper = Times(i+1);
        end

        KGt = IntKernelComp_TVHP( Times(1:i), option, upper, lower);

        Loglike = Loglike + sum(sum(KGt'.*Amc(Events(1:i),:)));

        if i>1
            tj = Times(1:i-1);
            cj = Events(1:i-1);


            acicj = model.A(cj,:,ci);
            gt = Kernel_TVHP(ti, tj, option.sigmaT, option.typeT);




            pijm = repmat(gt(:), [1,option.M])...
                .* acicj...
                .* repmat(kappat(:,i)', [i-1,1]);




            lambdai = lambdai + sum(pijm(:));



        end

        if lambdai>0
            Loglike = Loglike - log(lambdai);
        end

    end

    Loglike = Loglike + (Tstop-Tstart)*sum(model.mu(:));

end

Loglike = -Loglike;






