function Show_Intensity(Seqs, id, para, options)


History = [Seqs(id).Time; Seqs(id).Mark];
M = round(options.Tmax./options.dt);
time_stamp = 0:options.dt:(M-1)*options.dt;
lambda = zeros(size(para.A, 1), M);

for m = 1:M
    lambda(:, m) = Intensity_HP(time_stamp(m), History, para);
end
 
set(gcf,'Color','w');

figure
for u = 1:size(para.A, 1)
    subplot(1,size(para.A, 1),u)
    ind = find(History(2,:) == u);
    hold on
    stem(History(1, ind), ones(1, length(ind)), 'Color', 'k', 'LineWidth', 2);
    plot(time_stamp, lambda(u,:), 'Color', [0.8, 0, 0], 'LineWidth', 2);
    hold off
    ylabel('Intensity, \lambda(t)')
    xlabel(['Event-occurrence time (' num2str(length(ind)) ' events total)'])
end
    