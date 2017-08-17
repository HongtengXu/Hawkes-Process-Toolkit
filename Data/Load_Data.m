function [Seqs, Stats] = LoadData(file_path, file_format, ...
                                time_format, time_offset, time_scale)
                            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Load data and convert it to matlab format
%
% file_path: the path of data
% file_format: .csv or .txt file
% time_format: currently we support two formats. 1) real number; 2) real
% time, e.g., "year/month/day hour:minute:second" 3) ... you can define
% your format
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



Stats = struct('SeqID', [], ...
               'EventType', [], ...
               'FeatureType', [], ...
               'SeqLengthHist', [], ...
               'EventTypeHist', []);
          
if strcmp(file_format, 'txt')==0 && strcmp(file_format, 'csv')==0
    warning('The format might be unsupported.')
end

list = dir(sprintf('%s/*.%s', file_path, file_format));
if isempty(list)
    warning('None of supported files are found.')
end

Data = [];
num_file = length(list);
tic

if num_file>100
    parfor n = 1:num_file
        tic
        data = readtable(sprintf('%s/%s', file_path, list(n).name));
        Data = [Data; data];
        if mod(n, 100)==0 || n==num_file
            fprintf('File %d/%d, time=%.2fsec\n', n, num_file, toc);
        end
    end

else
    for n = 1:num_file

        data = readtable(sprintf('%s/%s', file_path, list(n).name));
        Data = [Data; data];
        if mod(n, 100)==0 || n==num_file
            fprintf('File %d/%d, time=%.2fsec\n', n, num_file, toc);
        end
    end
end


IDs = Data(:,1);
[userID, ~, newIDs] = unique(IDs);
Stats.SeqID = table2cell(userID);

Seqs = struct('Time', [], ...
              'Mark', [], ...
              'Start', [], ...
              'Stop', [], ...
              'Feature', []);
for i = 1:length(Stats.SeqID)
    Seqs(i).Start = 0;
end

Events = Data(:,3);
[eventType, ~, newEventType] = unique(Events);
Stats.EventType = table2cell(eventType);
Stats.EventTypeHist = zeros(size(eventType,1),1);

Time = table2cell(Data(:,2));

if size(Data,2)>3
    newFeatureType = [];
    Stats.FeatureType = cell(size(Data,2)-3,1);
    for i = 1:size(Data,2)-3
        Features = Data(:,i+3);
        [featureType, ~, newType] = unique(Features);
        Stats.FeatureType{i} = table2cell(featureType);
        newFeatureType = [newFeatureType, newType];        
    end    
end

for i = 1:size(Data,1)
    % sequence ID
    id = newIDs(i);
    
    % event ID
    event = newEventType(i);
    Stats.EventTypeHist(event) = Stats.EventTypeHist(event)+1;
    Seqs(id).Mark = [Seqs(id).Mark; event];
    
    % time stamp
    time = Time{i};
    switch time_format
        case 1 % real number
            t = (time-time_offset)/time_scale;

        case 2 % "year/month/day hour:minute:second, XXXX/XX/XX XX:XX:XX"
            Month = [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
            Month = cumsum(Month);
            year = str2double(time(1:4));
            month = str2double(time(6:7));
            day = str2double(time(9:10));
            hour = str2double(time(12:13));
            minute = str2double(time(15:16));
            second = str2double(time(18:19));
            
            t = (year-1)*365*24*60 + (Month(month)+day-1)*24*60 + ...
                hour*60 + minute + second/60;
            t = (t-time_offset)/time_scale;
            
        otherwise
            t = myTimeConvertFunc(time, time_offset, time_scale);
    end
    Seqs(id).Time = [Seqs(id).Time; t];
    
       
    if size(Data,2)>3
        Seqs(id).Feature = [Seqs(id).Feature; newFeatureType(i,:)];
    end  
    
    if mod(i, 10000)==0 || i == size(Data,1)
        fprintf('Events %d/%d, time=%.2fsec\n', i, size(Data,1), toc);
    end
end

L = zeros(length(Seqs),1);
if length(Seqs)>100
    parfor id = 1:length(Seqs)
        Seqs(id).Start = max([Seqs(id).Time(1)-eps, 0]);
        Seqs(id).Stop = Seqs(id).Time(end)+eps;
        L(id) = length(Seqs(id).Time);
    end
else
    for id = 1:length(Seqs)
        Seqs(id).Start = max([Seqs(id).Time(1)-eps, 0]);
        Seqs(id).Stop = Seqs(id).Time(end)+eps;
        L(id) = length(Seqs(id).Time);
    end
end

[counts, centers] = hist(L, min(L):max(L));
Stats.SeqLengthHist = [counts(:), centers(:)];

fprintf('Finish! Time=%.2fsec\n', toc);

end


function t = myTimeConvertFunc(time, time_offset, time_scale)
% users can define their own time convert function here!
end