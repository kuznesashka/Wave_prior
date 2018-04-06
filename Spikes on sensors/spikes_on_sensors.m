% -------------------------------------------------------
% Analysis for Jordy
% -------------------------------------------------------
%________________________________________________________
% Aleksandra Kuznetsova, kuznesashka@gmail.com
% Alexei Ossadtchi, ossadtchi@gmail.com

% 1. Import from bst Jordy data as Data

% 2. Automatic spike detection
f_low = 5;
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
% spike_mnl = Data.Events(1).times*Fs+1;
% clear ValMax IndMax
% for j = 1:length(spike_mnl)
%     spike = Ff(:,(spike_mnl(j)-10):(spike_mnl(j)+10));
%     %spike = Ff(:,(spike_mnl(j)));
%     [U,S,V] = svd(spike);
%     h = cumsum((diag(S)/norm(diag(S))).^2);
%     n = find(h>=0.98);
%     corr = MUSIC_scan(G2, U(:,1:n(1)));
%     
%     [ValMax(j), IndMax(j)] = max(corr);
%     j
% end

% figure
% hist(ValMax)

% ind_m = find(ValMax>=0.98);

% 3.2 for automatically detected spikes
clear ValMax IndMax
for j = 1:length(spike_ind)
    spike = Ff(:,((spike_ind(j)-10):(spike_ind(j)+10)));
    [U,S,V] = svd(spike);
    h = cumsum((diag(S)/norm(diag(S))).^2);
    n = find(h>=0.95);
    corr = MUSIC_scan(G2, U(:,1:n(1)));
    
    [ValMax(j), IndMax(j)] = max(corr);
    %loc{j} = RAP_MUSIC_scan(G2, U(:,1:n(1)), 0.97);
    j
end

figure
hist(ValMax)

ind_m = find(ValMax>=0.97);
spike_ind_m = spike_ind(ind_m);

IndMax_m = sort(IndMax(ind_m), 'ascend');
y = accumarray(IndMax_m(:),1);
y = y(y>0);
IndMax_m = unique(IndMax_m);

figure
h = scatter3(G3.GridLoc(1:50:length(G3.GridLoc(:,1)),1),...
    G3.GridLoc(1:50:length(G3.GridLoc(:,1)),2),...
    G3.GridLoc(1:50:length(G3.GridLoc(:,1)),3));
axis equal
grid off
axis off
view(360, 360)
hold on
scatter3(G3.GridLoc(IndMax_m,1),G3.GridLoc(IndMax_m,2),G3.GridLoc(IndMax_m,3), 100, y, 'filled')


% 4. Clustering
src_idx = IndMax(ind_m);

% stupid loop, there must be a different way
for i = 1:length(src_idx)
    for j = 1:length(src_idx)
        dist(i,j) = norm(G3.GridLoc(src_idx(i),:)-G3.GridLoc(src_idx(j),:));
    end
end

thr_dist = 0.01;
k = 1;
Nmin = 5;
fl = 1;
while fl == 1
    dst = sum(dist<thr_dist,2);
    [val, ind] = max(dst);
    if val > Nmin
        ind_nbh = find(dist(ind,:)<thr_dist);
        cluster{k} = src_idx(ind_nbh);
        spike_time{k} = spike_ind_m(ind_nbh);
        src_idx = src_idx(setdiff(1:size(dist,1), ind_nbh));
        spike_ind_m = spike_ind_m(setdiff(1:size(dist,1), ind_nbh));
        dist = dist(setdiff(1:size(dist,1), ind_nbh),setdiff(1:size(dist,1), ind_nbh));
        k = k + 1;
        fl = 1;
    else
        fl = 0;
    end
end

figure
h = scatter3(G3.GridLoc(1:50:length(G3.GridLoc(:,1)),1),...
    G3.GridLoc(1:50:length(G3.GridLoc(:,1)),2),...
    G3.GridLoc(1:50:length(G3.GridLoc(:,1)),3));
axis equal
grid off
axis off
view(360, 360)
hold on
for i = 1:length(cluster)
    scatter3(G3.GridLoc(cluster{i},1),G3.GridLoc(cluster{i},2), ...
        G3.GridLoc(cluster{i},3), 100, repmat(i, 1, length(cluster{i})), 'filled')
end


% 5. Phase analysis on sensors

% butterworth filter
[b,a] = butter(8,[25 35]/(Fs/2)); 
Ff = filtfilt(b,a,Data.F(32:99,:)')';
Ff_an = hilbert(Ff')';

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

for i = length(cluster)
    for j = 1:length(spike_time{i})
        figure
        plot_topo(X, Y, F_phase(:,spike_time{i}(j)))
        caxis([-10 10])
        colorbar
    end
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

