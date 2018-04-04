function [spike_ind_red] = SpikeDetect(Data, f_low, f_high, check)
% -------------------------------------------------------
% Automatical spike detection
% -------------------------------------------------------
% FORMAT:
%   SpikeDetect(Data)
% INPUTS:
%   Data -- Nch x Ntimes raw data 
%
% OUTPUT:
% 
% _______________________________________________________
% Aleksandra Kuznetsova, kuznesashka@gmail.com
% Alexei Ossadtchi, ossadtchi@gmail.com

Fs = 1/Data.Time(2)-Data.Time(1);
[b,a] = butter(3,[f_low f_high]/(Fs/2)); 
Ff = filtfilt(b,a,Data.F(32:99,:)')';

% 40 independent components
ncomp = 40;
[w,s] = runica(Ff, 'pca', ncomp);
W = w*s;
Q = pinv(W);
ica_ts = W*Ff;

% 2.1. Kurtosis (how outlier-prone the distribution is)
kurt = kurtosis(ica_ts');
[val_kurt, ind_kurt] = sort(kurt, 'descend');

if check == 1
    num_val = Fs*10;
    time = 1/Fs:1/Fs:(size(ica_ts,2)/Fs);
    ica_ts_sort = ica_ts(ind_kurt,:);
    spikes = int32(Data.Events(1).times*Fs+1);
    events = zeros(1, size(time, 2));
    events(spikes) = 1;

    scrolling_plot(time, ica_ts_sort, num_val, events, 2, val_kurt)
end

pos_ind = ind_kurt(1:20);

% 2.2. Spectral decomposition of components
for i = 1:length(pos_ind)
    S = fft(ica_ts(pos_ind(i),:));
    f = linspace(0, 250, length(S));
    Pyy = S.*conj(S)/length(S);
    
    power_low = sum(Pyy(f<=5));
    power_high = sum(Pyy((f>5)&(f<=125)));
    power_index(i) = (power_low/sum(f<=5))/(power_high/sum((f>5)&(f<=125)));
end

% figure
% stem(power_index)

% figure
% stem(log(power_index))

ind_power= find(log(power_index)<(mean(log(power_index))-std(log(power_index))));

% 3. Match filter
x = Ff(:,spikes(10)); % ideal spike
x = x/norm(x);
X_n = sqrt(sum(Ff.^2,1)); % norm of the each time sample
X = bsxfun(@rdivide, Ff, X_n);
match_filter = x'*X;
match_filter(match_filter<0.8) = 0;

ica_pos = ica_ts(pos_ind,:);
ica_pos_n = sqrt(sum(ica_pos.^2,2));
ica_pos = bsxfun(@rdivide, ica_pos, ica_pos_n);
match_filter = match_filter/norm(match_filter);
match_corr = sum(ica_pos.*match_filter,2);

% figure
% stem(match_corr)

ind_match= find(match_corr>(mean(match_corr)+std(match_corr)));

% Final decision
ind_spikecomp = intersect(ind_match, ind_power);

[d, ~] = max(abs(ica_ts(ind_spikecomp,:))); % maximum values
[val_d, ind_d] = sort(d); % threshold
thr = val_d(round(0.98*size(d,2))); 

% figure
% plot(1:size(d,2),d)
% hold on
% plot(1:size(d,2), repmat(thr,1, size(d,2)))

spike_ind = find(d>=thr);
% size(spike_ind)

diff = [30, spike_ind(2:length(spike_ind))-spike_ind(1:(length(spike_ind)-1))];
spike_ind_red = spike_ind(diff>20);
% size(spike_ind_red)

if check == 1
    print('Number of picked components ', num2str(ind_spikecomp))
    num_val = Fs*10;
    time = 1/Fs:1/Fs:(size(ica_ts,2)/Fs);
    spikes = int32(Data.Events(1).times*Fs+1);
    events = zeros(1, size(time, 2));
    events(spikes) = 1;
    events2 = zeros(1, size(time, 2));
    events2(spike_ind_red) = 1;

    scrolling_plot2(time, num_val, events, events2)
end

end
