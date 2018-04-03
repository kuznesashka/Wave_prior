% import from bst Jordy data as Data
Fs = 250; % sampling frequency
% butterworth filter
[b,a] = butter(8,[25 35]/(Fs/2)); 
Ff = filtfilt(b,a,Data.F(32:99,:)')';

Ff_an = hilbert(Ff')';
spike_ind = Data.Events.times/0.004+1;

k = 1;
for i = 32:99
    chan_loc(k,:) = mean(channel.Channel(i).Loc, 2)';
    k = k+1;
end


% 1. 
F_phase = unwrap(angle(Ff_an),2);

figure
scatter3(chan_loc(:,1),chan_loc(:,2),chan_loc(:,3), ...
    repmat(100,1,length(chan_loc(:,1))), F_phase(:,spike_ind(1)))

figure
imagesc(F_phase(:,(spike_ind(1)):(spike_ind(1)+50)))

[X,Y] = bst_project_2d(chan_loc(:,1), chan_loc(:,2), chan_loc(:,3), '2dlayout');
figure
plot_topo(X, Y, F_phase(:,spike_ind(1)))


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

