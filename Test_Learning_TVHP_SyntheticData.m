% demo: testing synthetic data

clear

load SynData_TVHP.mat


%type = 1;
%ShowTMHP(Seq{1,type}, mu, T, w, Period, Shift, MaxInfect, type)
SeqTrain = Seqs{2}(1:1000);
SeqTest = Seqs{2}(1001:2000);



option1.outer = 2;
option1.inner = 5;
option1.innerW = 5;
option1.thresW = 1e-5;
option1.step = 1e-5;
option1.D = 2;
option1.landmark = 0:5:50;
option1.M = length(option1.landmark);
option1.sigmaA = 5;
option1.typeA = 'gauss';
option1.sigmaT = 1;
option1.typeT = 'exp';
option1.thres = 1e-5;
option1.bias = 0.1;
option1.delta = 2;
option1.reg = 'L1';
option1.rho = 1;
option1.factor = 1.1;
option1.alphaS = 1;


NN = 5;
model = cell(NN,1);


    
for nn = 1:NN

    number = 200*nn;
    model{nn} = Learning_TVHP( SeqTrain(1:number), option1 );

end


save('Result_TVHP_Synthetic.mat', ...
    'model', 'option1');

%%
load Result_TVHP_Synthetic.mat


U = option1.D;
NN = 5;
A = zeros(U,U,T,Type);
Error = zeros(1, NN);
LogLike = zeros(1,NN);


for nn = 1:NN
       
    LogLike(nn) = Loglike_TVHP( SeqTest, model{nn}, option1 );
    Err = 0;
    for type = Type
        for t=1:T
            A(:,:,t,type) = Infectivity_TVHP(T, t, Period, Shift, MaxInfect, type);
        end

        basis = Kernel_TVHP( 1:T, option1.landmark, option1.sigmaA, option1.typeA );
        %figure
        for u=1:option1.D
            for v = 1:option1.D
                subplot(option1.D,option1.D,option1.D*(u-1)+v)
                tmp=A(u,v,:,type);
                Atmp = model{nn}.A(v,:,u)*basis;
                hold on
                plot(1:T, tmp(:), 'k-');
                plot(1:T, Atmp(:), 'r-');
                title(sprintf('a_{%d%d}(t)',u,v));
                legend('Real', 'Est')
                hold off
                Err = Err+norm(tmp(:)-Atmp(:));
            end
        end
        Error(nn) = Err/(option1.D^2);
    end
    
end



NUM = 200:200:1000;
figure
subplot(121)
hold on
plot(NUM, -LogLike, 'bs-');
axis tight
axis square
ylabel('LogLike');
xlabel('The number of training samples');
hold off

subplot(122)
hold on
plot(NUM, Error, 'bs-');
axis tight
axis square
ylabel('Relative Error');
xlabel('The number of training samples');
hold off