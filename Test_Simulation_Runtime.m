%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Test runtime of different simulation methods for Hawkes processes
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

D = 4; % the dimension of Hawkes processes
para1.mu = rand(D,1)/D;
para1.A = rand(D, D);
para1.A = 0.25 * para1.A./max(abs(eig(para1.A)));
para1.A = reshape(para1.A, [D, 1, D]);
para1.w = 1;
para2 = para1;
para2.kernel = 'exp';
para2.landmark = 0;
para3 = para2;

options.N = 100; % the number of sequences
options.Nmax = 1e4; % the maximum number of events per sequence
options.tstep = 0.2;
options.M = 50;


Tmax = 40:40:200;
Time = zeros(5,length(Tmax));

for i = 1:length(Tmax)
    options.Tmax = Tmax(i); % the maximum size of time window
    disp('Fast simulation of Hawkes processes with exponential kernel')
    tic
    Seqs1 = SimulationFast_Thinning_ExpHP(para1, options);
    Time(1,i) = toc;

    disp('Thinning-based simulation of Hawkes processes with arbitrary kernel')
    tic
    Seqs2 = Simulation_Thinning_HP(para2, options);
    Time(2,i) = toc;
    
    disp('Approximate simulation of Hawkes processes via branching process')
    tic
    options.GenerationNum = 3;
    Seqs3 = Simulation_Branch_HP(para3, options);
    Time(3,i) = toc;
    
    tic
    options.GenerationNum = 5;
    Seqs4 = Simulation_Branch_HP(para3, options);
    Time(4,i) = toc;
    
    tic
    options.GenerationNum = 10;
    Seqs5 = Simulation_Branch_HP(para3, options);
    Time(5,i) = toc;
end

figure
hold on
for i = 1:5
    plot(Tmax, log10(Time(i,:)), '-');
end
axis tight
xlabel('Length of time window');
ylabel('log Runtime (sec)');
legend('FastThin', 'Thin', 'Branch 3 Generations', ...
    'Branch 5 Generations', 'Branch 10 Generations')
