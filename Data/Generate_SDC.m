function Seqs_SDC = Generate_SDC(Seqs_Long, rate)

Seqs_SDC = Seqs_Long;

parfor n = 1:length(Seqs_Long)
    Start_new = Seqs_Long(n).Start + ...
        rand*(Seqs_Long(n).Stop - Seqs_Long(n).Start)*(1-rate);
    
    Stop_new = Start_new + rate*(Seqs_Long(n).Stop - Seqs_Long(n).Start);
    
    Seqs_SDC(n).Start = Start_new;
    Seqs_SDC(n).Stop = Stop_new;
    
    ind = find(Seqs_SDC(n).Time>=Start_new && Seqs_SDC(n).Time<=Stop_new);
    Seqs_SDC(n).Time = Seqs_SDC(n).Time(ind);
    Seqs_SDC(n).Mark = Seqs_SDC(n).Mark(ind);
    if ~isempty(Seqs_SDC(n).Feature)
        Seqs_SDC(n).Feature = Seqs_SDC(n).Feature(ind,:);
    end
end