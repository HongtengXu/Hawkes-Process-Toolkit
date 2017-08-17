function Show_Intensity_TVHP(Seq, mu, T, w, Period, Shift, MaxInfect, Type)

x = linspace(0,T,100);
k = length(mu);
%subplot(List)
set(gcf,'Color','w');
lambda = zeros(k, size(Seq.Time,2));
lambda(:,1)=mu;
for i = 2:size(Seq.Time,2)
    lambda(:,i) = Intensity_TVHP(Seq.Time(i), ...
        [Seq.Time;Seq.Mark], T, mu, w, Period, Shift, MaxInfect, Type);
end

figure

for i = 1:k
    subplot(1,k,i)
    
    %subplot(k,1,i)

    y = 0.05*max(lambda(i,:));
    e = find(Seq.Mark==i);
    for j = e
        line([Seq.Time(j) Seq.Time(j)],[0 y],'Color','k','LineWidth',2), hold on
    end    
    plot(Seq.Time,lambda(i,:),'Color',[.8 0 0],'LineWidth',2)
    axis([0 T 0 1.1*max(lambda(i,e))])
    ylabel(['Intensity, \lambda(t)'])
    xlabel(['Event-occurrence time (' num2str(length(e)) ' events total)'])
end