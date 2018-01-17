function [seqs, label] = SuperPosition(seqs1, seqs2, rate)


N = min([length(seqs1), length(seqs2)]);

seqs = struct('Time', [], ...
              'Mark', [], ...
              'Start', [], ...
              'Stop', [], ...
              'Feature', []);%'SourceLabel', []);
label = struct('source', []);
for n = 1:N
    if isempty(seqs1(n).Time)
        seqs(n).Time = seqs2(n).Time;
        seqs(n).Mark = seqs2(n).Mark;
        seqs(n).Start = seqs2(n).Start;
        seqs(n).Stop = seqs2(n).Stop;
        seqs(n).Feature = seqs2(n).Feature;
        label(n).source = 2*ones(size(seqs2(n).Time));
        %seqs(n).SourceLabel = 2*ones(size(seqs2(n).Time));
    else if isempty(seqs2(n).Time)
            seqs(n).Time = seqs1(n).Time;
            seqs(n).Mark = seqs1(n).Mark;
            seqs(n).Start = seqs1(n).Start;
            seqs(n).Stop = seqs1(n).Stop;
            seqs(n).Feature = seqs1(n).Feature;
            label(n).source = ones(size(seqs1(n).Time));
            %seqs(n).SourceLabel = ones(size(seqs2(n).Time));
        else
            Time1 = seqs1(n).Time;
            Mark1 = seqs1(n).Mark;
            Start1 = seqs1(n).Start;
            Stop1 = seqs1(n).Stop;

            Time2 = seqs2(n).Time;
            Mark2 = seqs2(n).Mark;
            Start2 = seqs2(n).Start;
            Stop2 = seqs2(n).Stop;

            dT = rate*(Stop1 - Start1);
            Start2_New = Start1 + dT;
            dT2 = Start2_New - Start2;

            Time2_New = Time2 + dT2;

            Time_New = [Time1, Time2_New];
            Mark_New = [Mark1, Mark2];
            tmp = [ones(size(Time1)), 2*ones(size(Time2))];
            
%             if ~isfield(seqs1(n), 'SourceLabel')
%                 Label1 = ones(1,length(Time1));
%             else
%                 Label1 = seqs1(n).SourceLabel;
%             end

%             if ~isfield(seqs2(n), 'SourceLabel')
%                 if ~isempty(Label1)
%                     value = max(Label1(:))+1;
%                 else
%                     value = 1;
%                 end
%                 Label2 = value*ones(1,length(Time2_New));
%             else
%                 if ~isempty(Label1)
%                     value = max(Label1(:))+1;
%                 else
%                     value = 1;
%                 end
%                 Label2 = value + ...
%                     (seqs2(n).SourceLabel - min(seqs2(n).SourceLabel));
%             end

            %Label = [Label1, Label2];
            Start_New = Start1;
            Stop_New = max([Stop1, Stop2 + dT2]);

            [Time_New, index] = sort(Time_New, 'ascend');

            Mark_New = Mark_New(index);
            label(n).source = tmp(index);
            seqs(n).Time = Time_New;
            seqs(n).Mark = Mark_New;
            seqs(n).Start = Start_New;
            seqs(n).Stop = Stop_New;
            %seqs(n).SourceLabel = Label(index);
        end
    end
end