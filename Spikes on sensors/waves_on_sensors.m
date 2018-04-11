% 5. Phase analysis on sensors
% for one spike

[b,a] = butter(8,[5 40]/(Fs/2)); 
Ff = filtfilt(b,a,Data.F(32:99,:)')';

spike_mnl = Data.Events(1).times*Fs+1;
ind = spike_mnl(5);

figure
plot(Ff(:, (ind-10):(ind+10))')

% phase for all channels
[b,a] = butter(8,[22 28]/(Fs/2)); 
Ff = filtfilt(b,a,Data.F(32:99,:)')';
Ff_an = hilbert(Ff')';

spike_phase = angle(Ff_an(:, ind:(ind+10)));

figure
stem(spike_phase(:,1))

% look at phases greater in absolute value than pi/2

spike_chan = find(abs(spike_phase(:,1)) > pi/2);
k = 1;
for i = 32:99
    chan_loc(k,:) = mean(channel.Channel(i).Loc, 2)';
    k = k+1;
end

[X,Y] = bst_project_2d(chan_loc(:,1), chan_loc(:,2), ...
    chan_loc(:,3),'2dlayout');

phase_latency = zeros(1, 68);
phase_latency(spike_chan) = (pi-abs(spike_phase(spike_chan,1)));

figure
plot_topo(X,Y,phase_latency')
caxis([min(phase_latency(spike_chan)) max(phase_latency)])

figure
scatter3(chan_loc(:,1),chan_loc(:,2),chan_loc(:,3), ...
    repmat(100,1,length(chan_loc(:,1))), phase_latency)

% crosscorrelation
for j = 1:68
r(:,j) = xcorr(Ff(1, ind:(ind+10)), ...
        Ff(j, ind:(ind+10)));
end
figure
imagesc(r)

figure
for i = 1:length(spike_chan)
    for j = 1:length(spike_chan)
        r(:,j) = xcorr(Ff(spike_chan(i), ind:(ind+10)), ...
        Ff(spike_chan(j), ind:(ind+10)));
    end
    subplot(5,3,i)
    imagesc(r)
    colorbar
end

figure
stem(-10:10,r)

% 6. For all spikes
% butterworth filter
[b,a] = butter(8,[22 28]/(Fs/2)); 
Ff = filtfilt(b,a,Data.F(32:99,:)')';
Ff_an = hilbert(Ff')';

k = 1;
for i = 32:99
    chan_loc(k,:) = mean(channel.Channel(i).Loc, 2)';
    k = k+1;
end

% 1. Only phase
F_phase = unwrap(angle(Ff_an),[],2);

% figure
% scatter3(chan_loc(:,1),chan_loc(:,2),chan_loc(:,3), ...
%     repmat(100,1,length(chan_loc(:,1))), F_phase(:,spike_mnl(1)))

% figure
% imagesc(F_phase(:,(spike_ind(1)):(spike_mnl(1)+50)))

[X,Y] = bst_project_2d(chan_loc(:,1), chan_loc(:,2), chan_loc(:,3),'2dlayout');

for i = length(cluster)
    for j = 1:length(cluster{i}(2,:))
        figure
        plot_topo(X, Y, F_phase(:,cluster{i}(2,j)))
        %caxis([-10 10])
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