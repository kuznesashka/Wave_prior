brainstorm
cd /home/sasha/WavePrior/git/waves

%% Initial parameters
PARAMS = struct();
PARAMS.protocol_name = 'Waves';
PARAMS.subject_name = 'Jordy_3';
PARAMS.study_name = 'jordy_spont03';
PARAMS.channels_idx = 32:99; 
PARAMS.sampling_rate = 250;
PARAMS.max_distance = 0.04; 
PARAMS.forward.name = 'Overlapping spheres_hr';
PARAMS.protocol_name = 'Spikes2';

PARAMS.cortex.low.comment = 'cortex_247306V';
PARAMS.visualize.flattening = false;
cortex_comment = PARAMS.cortex.low.comment;
cortex = from_bst_get_surface(cortex_comment, PARAMS);
faces = cortex.Faces;
vertices = cortex.Vertices;
PARAMS.wave.half_width = 0.005; % (m) half width of one wave
PARAMS.wave.duration = 0.02; % (s) time of wave duration 
PARAMS.wave.speeds = [0.001, 0.0025, 0.005, 0.01, 0.03 , 0.05, 0.07, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]; % (mm/ms = m/s)1

spike_time = {[11.488, 11.584],[12.116, 12.212], [41.940, 42.036], [42.160, 42.256], [42.640, 42.736]};
[meas, measfull]  = from_bst_get_measurements('jordy_spont03 (#1) band(1-40Hz)', PARAMS);
Ff = meas{1}.F(:,:);

% Forward model matrices
[G G3] = from_bst_get_gain_matrix(PARAMS.forward.name, PARAMS); 
Ns = size(G3,2)/3; % number of sources

ind = find(meas{1}.Time > 0.456 & meas{1}.Time < 0.5520); % averaged spike

figure
plot(Ff(:,ind)')

%% First point

Data = Ff(:,ind(1:5));
[U,S,V] = svd(Data);

MaxC = zeros(1,Ns);
SC = G'*U(:,1:7); 
for i=1:Ns
    MaxC(i) = sum(SC(i,:).^2); 
end;
[ValMax, IndMax] = max(MaxC); % the best fitted source, starting point

% force
IndMax = 198973;
figure
h = trimesh(faces,vertices(:,1),vertices(:,2),vertices(:,3), ones(1,size(vertices(:,1),1)));
set(h,'FaceAlpha',0.5);
axis equal
grid off
axis off
view(360, 360)
hold on
scatter3(vertices(IndMax,1),vertices(IndMax,2),vertices(IndMax,3), 100,'r', 'filled')

ind = find(meas{1}.Time > 0.468 & meas{1}.Time < 0.5240); % averaged spike
Data = Ff(:,ind);
figure
plot(Data')

%% Neighbourhood
DIST = zeros(1, size(vertices,1));
for d = 1:size(vertices,1)
    DIST(d) = norm(vertices(d,:)-vertices(IndMax,:));
end

ind_n = find(DIST<0.01);

surf_hr = cortex;
surf_lr = cortex;
vert_inds_lr = [];
vertices = cortex.Vertices;
faces = cortex.Faces;

minmsecort = repmat(0,1,size(vertices,1));
minmsecort(ind_n) = 1;
data_lr = minmsecort';

h = plot_brain_cmap2(surf_hr, surf_lr, vert_inds_lr, data_lr)
hold on
sphere_marker(vertices(IndMax, 1),vertices(IndMax, 2),vertices(IndMax, 3), 0.002, [0,0,0])

% sparse adjacency matrix 
PARAMS.VertConn = tess_vertconn(cortex.Vertices, cortex.Faces);
for npoint = 1:length(ind_n)
    PARAMS.seed_vertices = ind_n(npoint);
    sensor_waves{npoint} = wave_on_sensors_simple(cortex, PARAMS, G);
    npoint
end

%% localize the initial point
load sensor_waves_tot
load ind_n_tot

sensor_waves = sensor_waves_tot;
ind_n = ind_n_tot;

Nsens = size(Ff,1);
wavelength = PARAMS.wave.duration/(1/PARAMS.sampling_rate) + 1;
bestRsq = zeros(length(PARAMS.wave.speeds), length(ind_n));
bestshift = zeros(length(PARAMS.wave.speeds), length(ind_n));

for s = 1:length(PARAMS.wave.speeds)
    for npoint = 1:length(ind_n)
        rs = 1;
        clear DF SHIFT B STATS DataHat RSS Rsq wavesspeedlin waves wavesspeed
        waves = sensor_waves{npoint}; 
        Ndir = size(waves,1);
        wavesspeed = waves(:,s); 
        for d = 1:Ndir % for each direction
            tmp = wavesspeed{d}; % matrix (num sensors) by (num time points)
            wavesspeedlin(d,:) = tmp(:); % (num dir) by (num sensors)x(num time points)
        end

        R = size(Data,2)-wavelength; % time points minus the wave length
        range = 1:wavelength;
        for r = 1:R % sliding time window
            tmp = Data(:,range); % time interval with wave
            DataLin = tmp(:);
            TSS = sum((DataLin - mean(DataLin)).^2);
            lbd = lambda(wavesspeedlin',DataLin,'CV', 10, 'Alpha', 0.5);
            [B{rs},STATS{rs}] = lasso(wavesspeedlin',DataLin,'CV', 10, 'Alpha', 0.5, 'Lambda', lbd(70));
            DataHat{rs} = B{rs}'*wavesspeedlin +repmat(STATS{rs}.Intercept', 1, size(DataLin,1));
            SHIFT(rs,1) = r;
            diff = (DataHat{rs}'-repmat(DataLin,1,size(DataHat{rs},1)));
            RSS(rs,:) = sum((diff-repmat(mean(diff,1), size(diff,1), 1)).^2);
            Rsq(rs,:) = ones(1,size(DataHat{rs},1)) - RSS(rs,:)./TSS;
            range = range+1;
            rs = rs+1;
        end
      
        [bestRsq(s, npoint), bestshift(s, npoint)] = max(Rsq);
    end

    s
end

load bestRsq_tot
load ind_n_tot
load cortex_infl

ind_n = ind_n_tot;

s = 10;
%load cortex_infl
surf_hr_inf = cortex_infl2;
surf_hr = cortex;
surf_lr = cortex;
vert_inds_lr = [];
vertices = cortex.Vertices;
minmsecort = repmat(min(bestRsq(s,:)),1,size(vertices,1));
minmsecort(ind_n) = bestRsq(s,:);
data_lr = minmsecort';

% 208336 (747) s = 0.2
[val, idx] = max(bestRsq(s,:));
h = plot_brain_cmap2(surf_hr, surf_lr, vert_inds_lr, data_lr)
colorbar
hold on
scatter3(vertices(ind_n(idx),1),vertices(ind_n(idx),2),vertices(ind_n(idx),3), 50,'r', 'filled')

h = plot_brain_cmap2(surf_hr_inf, surf_lr, vert_inds_lr, data_lr)
colorbar
hold on
scatter3(cortex_infl.Vertices(ind_n(idx),1),cortex_infl.Vertices(ind_n(idx),2),cortex_infl.Vertices(ind_n(idx),3), 50,'r', 'filled')



%%
Nsens = size(Ff,1);
wavelength = PARAMS.wave.duration/(1/PARAMS.sampling_rate) + 1;
npoint = 1;
rs = 1;
clear DF MSE SPEED SHIFT B STATS DataHat RSS Rsq
for s = 1:length(PARAMS.wave.speeds);
    waves = sensor_waves{npoint}; 
    Ndir = size(waves,1);
    wavesspeed = waves(:,s); 
    for d = 1:Ndir % for each direction
        tmp = wavesspeed{d}; % matrix (num sensors) by (num time points)
        wavesspeedlin(d,:) = tmp(:); % (num dir) by (num sensors)x(num time points)
    end

    R = size(Data,2)-wavelength; % time points minus the wave length
    range = 1:wavelength;
    for r = 1:R % sliding time window
        tmp = Data(:,range); % time interval with wave
        DataLin = tmp(:);
        TSS = sum((DataLin - mean(DataLin)).^2);
        [B{rs},STATS{rs}] = lasso(wavesspeedlin(1:Ndir,:)',DataLin,'CV', 10, 'Alpha', 0.5);
        DataHat{rs} = B{rs}'*wavesspeedlin+repmat(STATS{rs}.Intercept', 1, size(DataLin,1));
        MSE(rs,:) = STATS{rs}.MSE;
        DF(rs,:) = STATS{rs}.DF;
        SPEED(rs) = PARAMS.wave.speeds(s);
        SHIFT(rs,1) = r;
        diff = (DataHat{rs}'-repmat(DataLin,1,size(DataHat{rs},1)));
        RSS(rs,:) = sum((diff-repmat(mean(diff,1), size(diff,1), 1)).^2);
        Rsq(rs,:) = ones(1,size(DataHat{rs},1)) - RSS(rs,:)./TSS;

        range = range+1;
        rs = rs+1;
    end
    
end
    

clear MaxRsq MaxRsqShift MinDF
range  = 1:R;
lambda = 1;
for s = 1:length(PARAMS.wave.speeds)
    [MaxRsq(s) MaxRsqShift(s)] = max(Rsq(range,lambda));
     MaxDF(s) = DF(range(MaxRsqShift(s)),lambda);
    range= range+R;
end;

figure
subplot(2,1,1)
plot(PARAMS.wave.speeds, MaxRsq, 'LineWidth', 1.5)
xlabel('Speed of propagation')
ylabel('Ratio of variance explained')
subplot(2,1,2)
plot(PARAMS.wave.speeds, MaxDF, 'LineWidth', 1.5)
xlabel('Speed of propagation')
ylabel('Degrees of freedom')

figure
imagesc(Rsq)

[valsort indsort] = sort(MaxRsq, 'descend');
bestspeed = PARAMS.wave.speeds(indsort(1));
bestshift = MaxRsqShift(indsort(1));
bestcoef = B{(indsort(1)-1)*R+bestshift}(:,30);
bestRsq = MaxRsq(indsort(1));

% picture
theta = 0:0.01:pi/4;
[val, ind] = sort(abs(bestcoef), 'descend');
bestcoef2 = bestcoef(ind);

figure
for i = 1:5
    if ind(i) ~=5
        rho = sin(2*theta).*cos(2*theta)*2*abs(bestcoef2(i));
        theta2 = theta+pi/8+pi/2*(ind(i)-3);
        if bestcoef2(i)>=0
            col = 'r';
        else col = 'b';
        end
        polar(theta2, rho, col)
        hold on
    else th = linspace(0,2*pi,50);
        r = abs(bestcoef2(i));
        if bestcoef2(i)>=0
            col = 'r';
        else col = 'b';
        end
  
    polar(th,r+zeros(size(th)), col)
    hold on
    end
end





