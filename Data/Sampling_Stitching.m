function Seqs_Stitched = Sampling_Stitching(Seqs_Short, para)


Seqs_Stitched = repmat(Seqs_Short, [1, para.ScaleUp]);

for i = 1:para.iteration
    for n = 1:length(Seqs_Stitched)
        weight = zeros(length(Seqs_Short), 2);
        for m = 1:length(Seqs_Short)
            if Seqs_Short(m).Stop < Seqs_Stitched(n).Start
                weight(m, 1) = Similarity(Seqs_Short(m).Stop, ...
                                          Seqs_Stitched(n).Start, ...
                                          para.sigma) * ...
                               Similarity(Seqs_Short(m).Feature, ...
                                          Seqs_Stitched(n).Feature, ...
                                          para.sigma);
            end
            
            if Seqs_Short(m).Start > Seqs_Stitched(n).Stop
                weight(m, 2) = Similarity(Seqs_Short(m).Start, ...
                                          Seqs_Stitched(n).Stop, ...
                                          para.sigma) * ...
                               Similarity(Seqs_Short(m).Feature, ...
                                          Seqs_Stitched(n).Feature, ...
                                          para.sigma);
            end
            
            if sum(weight(:,1))>0
                index = find(weight(:,1)>0);
                tmpw = cumsum(weight(index,1));
                value = rand * tmpw(end);
                for j = 1:length(tmpw)
                    if tmpw(j)>=value
                        break;
                    end
                end
                ind = index(j);
                
                Seqs_Stitched(n).Start = Seqs_Short(ind).Start;
                Seqs_Stitched(n).Time = [Seqs_Short(ind).Time; ...
                                        Seqs_Stitched(n).Time];
                Seqs_Stitched(n).Mark = [Seqs_Short(ind).Mark; ...
                                        Seqs_Stitched(n).Mark]; 
                Seqs_Stitched(n).Feature = [Seqs_Short(ind).Feature; ...
                                        Seqs_Stitched(n).Feature]; 
            end
            
            if sum(weight(:,2))>0
                index = find(weight(:,2)>0);
                tmpw = cumsum(weight(index,2));
                value = rand * tmpw(end);
                for j = 1:length(tmpw)
                    if tmpw(j)>=value
                        break;
                    end
                end
                ind = index(j);
                
                Seqs_Stitched(n).Stop = Seqs_Short(ind).Stop;
                Seqs_Stitched(n).Time = [Seqs_Stitched(n).Time; ...
                                        Seqs_Short(ind).Time];
                Seqs_Stitched(n).Mark = [Seqs_Stitched(n).Mark; ...
                                        Seqs_Short(ind).Mark]; 
                Seqs_Stitched(n).Feature = [Seqs_Stitched(n).Feature; ...
                                        Seqs_Short(ind).Feature]; 
            end
        end
        
    end
end

