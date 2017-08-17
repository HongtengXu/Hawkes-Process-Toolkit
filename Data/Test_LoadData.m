%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Test: load data and convert data format
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

time_format = 1;
time_offset = 0;
time_scale = 1;
[Seqs, Stats] = Load_Data('Linkedin_Data', 'txt', ...
                         time_format, time_offset, time_scale);
save('LinkedinData.mat', 'Seqs', 'Stats');

figure
subplot(121)
bar(Stats.SeqLengthHist(:,2),Stats.SeqLengthHist(:,1));
title('Histogram of sequence length');
subplot(122)
bar(1:length(Stats.EventTypeHist),Stats.EventTypeHist);
title('Histogram of event type');

time_format = 2;
time_offset = 2011*365*24*60;
time_scale = 60;
[Seqs, Stats] = Load_Data('IPTV_Data', 'txt', ...
                         time_format, time_offset, time_scale);
save('IPTVData.mat', 'Seqs', 'Stats');

figure
subplot(121)
bar(Stats.SeqLengthHist(:,2),Stats.SeqLengthHist(:,1));
title('Histogram of sequence length');
subplot(122)
bar(1:length(Stats.EventTypeHist),Stats.EventTypeHist);
title('Histogram of event type');