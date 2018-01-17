function model = Learning_WLS_HP(Seqs, model, alg, flag)

D = size(model.A, 1);
L = length(model.landmark);
C = length(Seqs);
XX = [];
YY = [];

for c = 1:C
    
    Time = Seqs(c).Time;
    if ~isempty(Time)
    Event = Seqs(c).Mark;
    Tstart = Seqs(c).Start;

%     if isempty(alg.Tmax)
%         Tstop = Seqs(c).Stop;
%     else
%         Tstop = alg.Tmax;
%         indt = Time < alg.Tmax;
%         Time = Time(indt);
%         Event = Event(indt);
%     end

    Y = zeros(length(Time), 1);
    
    if strcmp(flag, 'multisource')
        Xmu = zeros(length(Time), D*C);
    else
        Xmu = zeros(length(Time), D);
    end
    XA = zeros(length(Time), D*L*D);
    
    for i = 1:length(Time)
        di = Event(i);
        if strcmp(flag, 'multisource')
            Xmu(i, D*(c-1)+di) = 1;
        else
            Xmu(i, di)= 1;
        end
        
        index = find(Event(1:i)==di);
        Y(i) = length(index);
        
        if i>1
            dt = Time(i) - Time(1:i-1);
            G = Kernel_Integration(dt, model);
            
            for j = 1:i-1
                dj = Event(j);
                Gij = G(j,:);
                XA(i, D*L*(di-1)+(dj:L:D*L)) = ...
                    XA(i, D*L*(di-1)+(dj:L:D*L)) + Gij';
            end
        end
    end
    
    Ytmp = alg.w.*Y./sqrt((Time(:) - Tstart));
    Xtmp = alg.w.*[Xmu, XA]./repmat(sqrt(Time(:) - Tstart), ...
        [1, size(Xmu,2)+size(XA,2)]);
    YY = [YY; Ytmp];
    XX = [XX; Xtmp];
    end
end

res0 = rand(size(XX,2),1);
options = optimoptions('fmincon', 'Algorithm', 'interior-point', ...
    'Display', 'iter', ...
    'MaxFunctionEvaluations', 30000);
res = fmincon(@(x)norm(XX*x - YY)^2, res0, [], [], [], [], ...
    zeros(size(XX,2),1), [], [], options);
%res = XX\YY;
if strcmp(flag, 'multisource')
    model.mu = reshape(res(1:C*D), [D, C]);
    model.A = reshape(res(1+C*D:end), [D, L, D]);
else
    model.mu = reshape(res(1:D), [D, 1]);
    model.A = reshape(res(1+D:end), [D, L, D]);
end