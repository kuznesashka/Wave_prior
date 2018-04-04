% -------------------------------------------------------
% Analysis for Jordy
% -------------------------------------------------------
%________________________________________________________
% Aleksandra Kuznetsova, kuznesashka@gmail.com
% Alexei Ossadtchi, ossadtchi@gmail.com

% 1. Import from bst Jordy data as Data
% 2. Automatic spike detection

f_low = 1;
f_high = 45;
spike_ind = SpikeDetect(Data, f_low, f_high, 1);

Fs = 1/Data.Time(2)-Data.Time(1);
[b,a] = butter(3,[f_low f_high]/(Fs/2)); 
Ff = filtfilt(b,a,Data.F(32:99,:)')';

% 3. Spikes source localization
% Import from bst G3 forward matrix
G3.Gain = G3.Gain(32:99,:);
Nsites = size(G3.Gain,2)/3;
Ns = size(G3.Gain,1);

range3 = 1:3;
range2 = 1:2;
G2 = zeros(size(G3.Gain,1),2*Nsites);
for i = 1:Nsites
    g = G3.Gain(:,range3);
    [u s v] = svd(g, 'econ');
    G2(:,range2) = u(:,1:2);
    range3 = range3+3;
    range2 = range2+2;
end

% 3.1 for manually detected spikes
spike_mnl = Data.Events(1).times*Fs+1;
clear ValMax IndMax
for j = 1:length(spike_mnl)
    spike = Ff(:,(spike_mnl(j)-10):(spike_mnl(j)+10));
    [U,S,V] = svd(spike);
    h = cumsum((diag(S)/norm(diag(S))).^2);
    n = find(h>=0.95);
    corr = MUSIC_scan(G2, U(:,1:n(1)));
    
    [ValMax(j), IndMax(j)] = max(corr);
    j
end

figure
hist(ValMax)

ind_m = find(ValMax>=0.97);

figure
h = scatter3(G3.GridLoc(1:20:length(G3.GridLoc(:,1)),1),...
    G3.GridLoc(1:20:length(G3.GridLoc(:,1)),2),...
    G3.GridLoc(1:20:length(G3.GridLoc(:,1)),3));
axis equal
grid off
axis off
view(360, 360)
hold on
scatter3(G3.GridLoc(IndMax(ind_m),1),G3.GridLoc(IndMax(ind_m),2),...
    G3.GridLoc(IndMax(ind_m),3), 100,'r', 'filled')


% 3.2 for automatically detected spikes
clear ValMax IndMax
for j = 1:length(spike_ind)
    spike = Ff(:,(spike_ind(j)-10):(spike_ind(j)+10));
    [U,S,V] = svd(spike);
    h = cumsum((diag(S)/norm(diag(S))).^2);
    n = find(h>=0.95);
    corr = MUSIC_scan(G2, U(:,1:n(1)));
    
    [ValMax(j), IndMax(j)] = max(corr);
    j
end

figure
hist(ValMax)

ind_m = find(ValMax>=0.97);

figure
h = scatter3(G3.GridLoc(1:20:length(G3.GridLoc(:,1)),1),...
    G3.GridLoc(1:20:length(G3.GridLoc(:,1)),2),...
    G3.GridLoc(1:20:length(G3.GridLoc(:,1)),3));
axis equal
grid off
axis off
view(360, 360)
hold on
scatter3(G3.GridLoc(IndMax(ind_m),1),G3.GridLoc(IndMax(ind_m),2),G3.GridLoc(IndMax(ind_m),3), 100,'r', 'filled')


% 4. Clustering
src_idx = IndMax(ind_m);

% stupid loop, there must be a different way
for i = 1:length(src_idx)
    for j = 1:length(src_idx)
        dist(i,j) = norm(G3.GridLoc(src_idx(i))-G3.GridLoc(src_idx(j)));
    end
end

thr_dist = 0.02;
dst = sum(dist>thr_dist);
[ind, val] = max(dst);
cluster{1} = 




% 5. Phase analysis on sensors

% butterworth filter
[b,a] = butter(8,[25 35]/(Fs/2)); 
Ff = filtfilt(b,a,Data.F(32:99,:)')';
Ff_an = hilbert(Ff')';

spike_mnl = Data.Events(1).times/0.004+1;

k = 1;
for i = 32:99
    chan_loc(k,:) = mean(channel.Channel(i).Loc, 2)';
    k = k+1;
end

% 1. Only phase
F_phase = unwrap(angle(Ff_an),2);

% figure
% scatter3(chan_loc(:,1),chan_loc(:,2),chan_loc(:,3), ...
%     repmat(100,1,length(chan_loc(:,1))), F_phase(:,spike_mnl(1)))

% figure
% imagesc(F_phase(:,(spike_ind(1)):(spike_mnl(1)+50)))

[X,Y] = bst_project_2d(chan_loc(:,1), chan_loc(:,2), chan_loc(:,3),'2dlayout');

for i = 1:length(spike_mnl)
    figure
    plot_topo(X, Y, F_phase(:,spike_mnl(i)))
    caxis([0 20])
    colorbar
end


% 2.
F_phase = angle(Ff_an(:,spike_ind(1):(spike_ind(1)+50)));
F_phase = (F_phase>0).*2*pi - F_phase; 

figure
imagesc(F_phase)
colorbar

figure
scatter3(chan_loc(:,1),chan_loc(:,2),chan_loc(:,3), ...
    repmat(100,1,length(chan_loc(:,1))), F_phase(:,1))

figure
plot_topo(X, Y, F_phase(:,1))

