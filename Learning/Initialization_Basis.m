function model = Initialization_Basis(Seqs, baseType, bandwidth, landmark)

D = zeros(length(Seqs),1);
for i = 1:length(Seqs)            
    D(i) = max(Seqs(i).Mark);
end
D = max(D);

switch nargin
    case 1
        sigma = zeros(D);
        Tmax = zeros(D);

        est = cell(D);
        id = randperm(length(Seqs));
        for n = 1:min([length(id), 10])
            for i = 2:length(Seqs(id(n)).Time)
                ti = Seqs(id(n)).Time(i);
                di = Seqs(id(n)).Mark(i);
                for j = 1:i-1
                    tj = Seqs(id(n)).Time(j);
                    dj = Seqs(id(n)).Mark(j);
                    est{di, dj} = [est{di, dj}, ti - tj];
                end
            end
        end
        
        for di = 1:D
            for dj = 1:D
                sigma(di, dj) = ((4*std(est{di, dj})^5)/(3*length(est{di, dj})))^0.2; 
                Tmax(di, dj) = mean(est{di, dj});
            end
        end
        Tmax = min(Tmax(:))/2;
        
        
        model.kernel = 'gauss';
        model.w = min(sigma(:))/2;
        model.landmark = model.w*(0:ceil(Tmax/model.w));
        
    case 2
        model.kernel = baseType;
        
        sigma = zeros(D);
        Tmax = zeros(D);

        est = cell(D);
        id = randperm(length(Seqs));
        for n = 1:min([length(id), 10])
            for i = 2:length(Seqs(id(n)).Time)
                ti = Seqs(id(n)).Time(i);
                di = Seqs(id(n)).Mark(i);
                for j = 1:i-1
                    tj = Seqs(id(n)).Time(j);
                    dj = Seqs(id(n)).Mark(j);
                    est{di, dj} = [est{di, dj}, ti - tj];
                end
            end
        end
        
        for di = 1:D
            for dj = 1:D
                sigma(di, dj) = ((4*std(est{di, dj})^5)/(3*length(est{di, dj})))^0.2; 
                Tmax(di, dj) = mean(est{di, dj});
            end
        end
        Tmax = min(Tmax(:))/2;        
        model.w = min(sigma(:))/2;
        model.landmark = model.w*(0:ceil(Tmax/model.w));
        
    case 3
        model.kernel = baseType;
        model.w = bandwidth;
        sigma = zeros(D);
        Tmax = zeros(D);

        est = cell(D);
        id = randperm(length(Seqs));
        for n = 1:min([length(id), 10])
            for i = 2:length(Seqs(id(n)).Time)
                ti = Seqs(id(n)).Time(i);
                di = Seqs(id(n)).Mark(i);
                for j = 1:i-1
                    tj = Seqs(id(n)).Time(j);
                    dj = Seqs(id(n)).Mark(j);
                    est{di, dj} = [est{di, dj}, ti - tj];
                end
            end
        end
        
        for di = 1:D
            for dj = 1:D
                sigma(di, dj) = ((4*std(est{di, dj})^5)/(3*length(est{di, dj})))^0.2; 
                Tmax(di, dj) = mean(est{di, dj});
            end
        end
        Tmax = min(Tmax(:))/2;
        
        model.landmark = model.w*(0:ceil(Tmax/model.w));
        
    otherwise
        model.kernel = baseType;
        model.w = bandwidth;
        model.landmark = landmark;
end



model.A = rand(D, length(model.landmark), D)./(D^2 * length(model.landmark));
model.mu = rand(D, 1)./D;