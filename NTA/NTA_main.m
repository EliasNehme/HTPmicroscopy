%% NTA: track linking, msd calculation, msd fitting

% start with a clean slate
close all; clear all; clc;

%% localizations and settings

% whether to show the tracking analysis for low/high density
% str_density = 'low_density';
str_density = 'high_density';

% load the deepstorm3d csv results
locs = importdata(char(strcat('localizations_deepstorm3d_',str_density,'.csv')));
fxyzi = locs.data;

% load all locations and confidences
f = fxyzi(:,1);
xyz = fxyzi(:,2:4)./1e3; % nm -> um
conf = fxyzi(:,5);

% thresholded localizations
T = 40.0;
ind = conf >= T & f <= 100.0;
f = f(ind);
xyz = xyz(ind,:);
conf = conf(ind);
num_frames = max(f); % should be 100 frames

% thresholds for track linking between successive frames in xy and in z
thresh_xy = 20;
thresh_z = 30;

% track length threshold (50% of total time lapse) for msd calculation
len_thresh = 50;

% experiment time settings
dt = 2.0; % time in between successive frames in [sec]
t_exp = 0.4; % exposure time in [sec]

% fraction of msd pts to fit with a linear line
frac = 0.02;

%% Sequential track linking

% go over frames and regenerate the images and look at the 3D positions
tracks = {};
active_tracks = [];
ntracks = 0;
for i=1:num_frames
    
    % current number of emitters and positions
    xyzi = xyz(f==i,:);
    confi = conf(f==i);
    num_emitters = size(xyzi,1);
    
    % if this is the first frame start tracks
    if i==1
        
        % start tracks
        for j=1:num_emitters
            ntracks = ntracks + 1;
            tracks{ntracks} = [i, xyzi(ntracks,:), conf(ntracks)]; % [frm #, x, y, z, conf]
            active_tracks(ntracks) = 1;
        end
    
    % otherwise match points and update tracks accordingly
    else
        
        % previous xyz is the last point in the tracks that are active
        xyz_prev = [];
        conf_prev = [];
        track_ind = [];
        for t=1:ntracks
            if active_tracks(t)
                track = tracks{t};
                xyz_prev = [xyz_prev; track(end,2:4)];
                conf_prev = [conf_prev, track(end,5)];
                track_ind = [track_ind, t];
            end
        end
        
        % calculate the distance matrix in 3D: D(i,j) = d(prev_i, new_j)
        D = pdist2(xyz_prev, xyzi, 'euclidean');
        Dxy = pdist2(xyz_prev(:,1:2), xyzi(:,1:2), 'euclidean');
        Dz = pdist2(xyz_prev(:,3), xyzi(:,3), 'euclidean');
        
        % modify distance in 3D to be infinite for unreasonable linkages
        D(Dxy(:) > thresh_xy) = Inf;
        D(Dz(:) > thresh_z) = Inf;

        % apply the hungarian algorithm vectorized implementation
        addpath('Hungarian Algorithm');
        matches = munkres(D);
        matches(matches==0) = NaN;
        
        % matched previous points and new points
        matched_prev = find(~isnan(matches));
        matched_curr = matches(matched_prev);
        
        % for matched previous tracks elongate them
        matched_tracks = track_ind(matched_prev);
        for j=1:length(matched_prev)
            matched_track = track_ind(matched_prev(j));
            curr_track = tracks{matched_track};
            match_ind_new = matched_curr(j);
            new_pt = xyzi(match_ind_new,:);
            new_conf = conf(match_ind_new);
            tracks{matched_track} = [curr_track; [i, new_pt, new_conf]];
        end
            
        % for new pts start new tracks
        new_pts = setdiff(1:num_emitters, matched_curr);
        for n=1:length(new_pts)
            ntracks = ntracks + 1;
            tracks{ntracks} = [i, xyzi(new_pts(n),:), conf(new_pts(n))]; % [frm #, x, y, z, conf]
            active_tracks(ntracks) = 1;
        end
    end

    % plot the current and compare to previous 3D positions if applicable
    figure(103);
    for t=1:ntracks
        track_i = tracks{t};
        plot3(track_i(:,2),track_i(:,3),track_i(:,4),'LineWidth',2);
        hold on;
    end
    hold off;
    title(['# of continuous tracks is ' num2str(ntracks) ' , frame # ' num2str(i)]);
    xlim([0,round(max(xyz(:,1)))+20]);
    ylim([0,round(max(xyz(:,2)))+20]);
    zlim([0,100]);
    set(gca,'Ydir','reverse');
    view(-10,60);    
    drawnow;
    % pause;
end

% close accumulating figure
close 103;

%% sorting by length and filtering

% sort by length
L = [];
for i=1:ntracks
    ti = tracks{i};
    L(i) = size(ti,1);
end
[Ls,I] = sort(L,'descend');

% reorder tracks: tall -> short
tracks_sort = {};
for i=1:ntracks
    tracks_sort{i} = double(tracks{I(i)});
end

% examine all tracks visually
hf = figure(104);
hf.Position = [800, 300, 700, 700];
cm = jet(100);
for i=1:ntracks
    track_i = tracks_sort{i};
    surface([track_i(:,2),track_i(:,2)],[track_i(:,3),track_i(:,3)],...
        [track_i(:,4),track_i(:,4)],[track_i(:,1),track_i(:,1)],...
        'facecol','no','edgecol','interp','linew',2);
    hold on;
    title(['Track # ' num2str(i) ' , track length is ' num2str(Ls(i)) ' frames']);
    xlim([0,round(max(xyz(:,1)))+20]);
    ylim([0,round(max(xyz(:,2)))+20]);
    zlim([0,100]);
    set(gca,'Ydir','reverse');
    if i==0
        view(0,90);
    else
        [caz, cel] = view;
        view(caz, cel);
    end
    colormap(cm);
    drawnow;
    if Ls(i)<7
        break;
    end
    % pause;
end

% close accumulating figure
close 104;

% save linked tracks after sorting
tracks = tracks_sort;
save(strcat('tracks_TP_',str_density,'.mat'), 'tracks');

%% threshold track length to exclude short tracks from the analysis

% record track lengths
num_tracks = length(tracks);
track_lens = zeros(num_tracks,1);
for i=1:num_tracks
    track_lens(i) = length(tracks{i});
end

% threshold short tracks
tracks2 = {}; k = 1;
for i=1:length(tracks)
    if size(tracks{i},1) >= len_thresh
        tracks2{k} = tracks{i};
        k = k + 1;
    end
end
tracks = tracks2;
ntracks = length(tracks);

% number of tracks above threshold
ind_aboveT = track_lens>=len_thresh;
num_tracksT = sum(ind_aboveT);
avg_len = mean(track_lens(ind_aboveT));
figure();
hist(track_lens(ind_aboveT),5);
hx = xlabel('Track Length');
hy = ylabel('Counts');
ht = title(['For ' num2str(num_tracksT) '/' num2str(length(track_lens)) ' tracks, Average Length is ' num2str(avg_len)]);
set([hx,hy,ht],'interpreter','latex','FontSize',15);

%% plot considered tracks with tips at starting point

% colormap for tracks and number of frames
cm_time = jet(num_frames);
cm_im = gray(num_frames);
cm = [cm_time; cm_im];

% examine all tracks visually
hf = figure(104);
hf.Position = [650, 300, 850, 410];
% hf.Position = [800, 300, 700, 700];
for i=1:ntracks
    
    % if below length threshold stop plotting
    if Ls(i)<len_thresh
        break;
    end
    
    % current track
    track_i = double(tracks{i});
    t = track_i(:,1)./num_frames*0.35 + 0.1;
    x = smooth(track_i(:,2),3);
    y = smooth(track_i(:,3),3);
    z = smooth(track_i(:,4),3);
    
    % plot current track in 3D
    surface([x,x],[y,y],[z,z],[t,t],'facecol','no','edgecol','interp','linew',2);
    hold on;
    
    % plot track tip at start
    scatter3(x(1),y(1),z(1),50,'black','filled');
    hold on;
    
    % add labels, title, axis limits and colormap
    xlabel(['x (' char(956) 'm)']);
    ylabel(['y (' char(956) 'm)']);
    zlabel(['z (' char(956) 'm)']);
    % title(['Track # ' num2str(i) ' , track length is ' num2str(L(i)) ' frames']);
    xlim([0, round(max(xyz(:,1)))+20]);
    ylim([0, round(max(xyz(:,2)))+20]);
    zlim([0, 100]);
    set(gca,'Ydir','reverse');
    grid on;
    colormap(cm);
    caxis([0, 1]);
    set(gca,'FontSize',15);
    daspect([1,1,1]);
    
    % set view
    view(-85,15);
    if i==0
        view(0,90);
    else
        [caz, cel] = view;
        view(caz, cel);
    end
    
    % pause for checking
    drawnow;
    % pause;
end

% more spaced ticks
set(gca,'XTick',[0,100,200,300],'XTickLabel',[0,100,200,300]);
set(gca,'YTick',[0,100,200,300],'YTickLabel',[0,100,200,300]);
set(gca,'ZTick',[0, 50, 100],'ZTickLabel',[0, 50, 100]);

% align xaxis and yaxis labels
shifts_x = [-50, +30, -30];
shifts_y = [-10, 0, -10];
FixLabelsPosition(gca,-85,15,shifts_x,shifts_y);

%% Calculate ensemble MSDs per axis

% time lags in seconds
tau_vec = (1:num_frames)*(dt+t_exp);

% make all tracks length num_frames using nans
tracks_nans ={};
for i=1:ntracks
    
    % current track
    track_i = tracks{i};
    ti = track_i(:,1);
    
    % make track with length 100 using nans
    track_i_nans = [];
    for f=1:num_frames
        indf = find(ti==f);
        if ~isempty(indf)
            track_i_nans = [track_i_nans;track_i(indf,:)];
        else
            track_i_nans = [track_i_nans;nan(1,5)];
        end
    end
    tracks_nans{i} = track_i_nans;
end

% for each tau calculate ensemble MSD
ensemble_msd = zeros(num_frames,11);
for tau=1:num_frames
    
    % for all possible pairs with lag tau calculate SD
    num_pairs_ensemble = 0;
    square_dist3d_ensemble = [];
    square_dist2d_ensemble = [];
    square_distx_ensemble = [];
    square_disty_ensemble = [];
    square_distz_ensemble = [];
    
    % for every track take all point pairs that are available for MSD(tau)
    for i=1:ntracks
        
        % geth ith track with length num_frames
        track_i_nans = tracks_nans{i};
    
        % for every tau get possible pair 
        for j=1:(num_frames-tau)

            % if the jth point and the (j-tau)th point are sampled
            % calculate the square distance
            if ~isnan(track_i_nans(j,1)) && ~isnan(track_i_nans(j+tau,1))
                num_pairs_ensemble = num_pairs_ensemble + 1;
                square_dist3d_ensemble(num_pairs_ensemble) = sum((track_i_nans(j,2:4) - track_i_nans(j+tau,2:4)).^2);
                square_dist2d_ensemble(num_pairs_ensemble) = sum((track_i_nans(j,2:3) - track_i_nans(j+tau,2:3)).^2);
                square_distx_ensemble(num_pairs_ensemble) = sum((track_i_nans(j,2) - track_i_nans(j+tau,2)).^2);
                square_disty_ensemble(num_pairs_ensemble) = sum((track_i_nans(j,3) - track_i_nans(j+tau,3)).^2);
                square_distz_ensemble(num_pairs_ensemble) = sum((track_i_nans(j,4) - track_i_nans(j+tau,4)).^2);
            end
        end
    end
    
    % save ensemble MSD
    ensemble_msd(tau,1) = mean(square_dist3d_ensemble);
    ensemble_msd(tau,2) = std(square_dist3d_ensemble);
    ensemble_msd(tau,3) = num_pairs_ensemble; 
    ensemble_msd(tau,4) = mean(square_dist2d_ensemble);
    ensemble_msd(tau,5) = std(square_dist2d_ensemble);
    ensemble_msd(tau,6) = mean(square_distx_ensemble);
    ensemble_msd(tau,7) = std(square_distx_ensemble);
    ensemble_msd(tau,8) = mean(square_disty_ensemble);
    ensemble_msd(tau,9) = std(square_disty_ensemble);
    ensemble_msd(tau,10) = mean(square_distz_ensemble);
    ensemble_msd(tau,11) = std(square_distz_ensemble);
    
    % report progress
    disp(['Finished calculating ensemble MSD for tau ' num2str(tau) ' /' num2str(num_frames)]);
end

% ensemble +- std
ensemble_msd_xyz = ensemble_msd(:,1)';
ensemble_msd_plus_std_xyz = ensemble_msd(:,1)' + ensemble_msd(:,2)';
ensemble_msd_minus_std_xyz = ensemble_msd(:,1)' - ensemble_msd(:,2)';
ensemble_msd_xy = ensemble_msd(:,4)';
ensemble_msd_plus_std_xy = ensemble_msd(:,4)' + ensemble_msd(:,5)';
ensemble_msd_minus_std_xy = ensemble_msd(:,4)' - ensemble_msd(:,5)';
ensemble_msd_x = ensemble_msd(:,6)';
ensemble_msd_plus_std_x = ensemble_msd(:,6)' + ensemble_msd(:,7)';
ensemble_msd_minus_std_x = ensemble_msd(:,6)' - ensemble_msd(:,7)';
ensemble_msd_y = ensemble_msd(:,8)';
ensemble_msd_plus_std_y = ensemble_msd(:,8)' + ensemble_msd(:,9)';
ensemble_msd_minus_std_y = ensemble_msd(:,8)' - ensemble_msd(:,9)';
ensemble_msd_z = ensemble_msd(:,10)';
ensemble_msd_plus_std_z = ensemble_msd(:,10)' + ensemble_msd(:,11)';
ensemble_msd_minus_std_z = ensemble_msd(:,10)' - ensemble_msd(:,11)';

% saving everything
save(strcat('ensemble_msds_TP_',str_density,'.mat'),'tau_vec','track_lens',...
    'ensemble_msd_xyz', 'ensemble_msd_plus_std_xyz', 'ensemble_msd_minus_std_xyz',...
    'ensemble_msd_xy', 'ensemble_msd_plus_std_xy', 'ensemble_msd_minus_std_xy',...
    'ensemble_msd_x', 'ensemble_msd_plus_std_x', 'ensemble_msd_minus_std_x',...
    'ensemble_msd_y', 'ensemble_msd_plus_std_y', 'ensemble_msd_minus_std_y',...
    'ensemble_msd_z', 'ensemble_msd_plus_std_z', 'ensemble_msd_minus_std_z');

%% Fit ensemble MSD

% fraction of msd points to fit
num_pts = length(tau_vec);
tau_fit = tau_vec(1:num_pts*frac);
tau_msd_0 = [0, tau_fit];

% fit XYZ
ensemble_fit_xyz = ensemble_msd_xyz(1:num_pts*frac);
pxyz = polyfit(tau_fit, ensemble_fit_xyz, 1);
ensemble_msd_0_xyz = polyval(pxyz, tau_msd_0);
D_xyz = pxyz(1)/2;
std_xyz = 1/sqrt(2)*sqrt(pxyz(2) + pxyz(1)*t_exp/3);

% only XY
ensemble_fit_xy = ensemble_msd_xy(1:num_pts*frac);
pxy = polyfit(tau_fit, ensemble_fit_xy, 1);
ensemble_msd_0_xy = polyval(pxy, tau_msd_0);
D_xy = pxy(1)/2;
std_xy = 1/sqrt(2)*sqrt(pxy(2) + pxy(1)*t_exp/3);

% only X
ensemble_fit_x = ensemble_msd_x(1:num_pts*frac);
px = polyfit(tau_fit, ensemble_fit_x, 1);
ensemble_msd_0_x = polyval(px, tau_msd_0);
D_x = px(1)/2;
std_x = 1/sqrt(2)*sqrt(px(2) + px(1)*t_exp/3);

% only Y
ensemble_fit_y = ensemble_msd_y(1:num_pts*frac);
py = polyfit(tau_fit, ensemble_fit_y, 1);
ensemble_msd_0_y = polyval(py, tau_msd_0);
D_y = py(1)/2;
std_y = 1/sqrt(2)*sqrt(py(2) + py(1)*t_exp/3);

% only Z
ensemble_fit_z = ensemble_msd_z(1:num_pts*frac);
pz = polyfit(tau_fit, ensemble_fit_z, 1);
ensemble_msd_0_z = polyval(pz, tau_msd_0);
D_z = pz(1)/2;
std_z = 1/sqrt(2)*sqrt(pz(2) + pz(1)*t_exp/3);

% ensemble MSD xyz
hf = figure();
hf.Position = [100,100, 1200, 800];
subplot(2,2,1);
pxx = [tau_vec(1:end-1), fliplr(tau_vec(1:end-1))];
pyy = [ensemble_msd_plus_std_xyz(1:end-1), fliplr(ensemble_msd_minus_std_xyz(1:end-1))];
fill(pxx, pyy, 1, 'FaceColor', [0.5,0.5,0.5],'EdgeColor', 'none', 'FaceAlpha',0.2);
hold on;
plot(tau_vec, ensemble_msd_xyz,'b','Linewidth',2);
hold on;
plot(tau_fit, ensemble_fit_xyz,'or','Linewidth',2);
hold on;
plot(tau_msd_0, ensemble_msd_0_xyz,'--k','Linewidth',2);
hold off;  
grid on;
xlim([0, 30]);
hx = xlabel('$\tau \ \left[Sec\right]$');
hy = ylabel('MSD $\left[\frac{\mu m^2}{sec}\right]$');
ht = title(['Ensemble MSD XYZ: D=' num2str(D_xyz,2) ', $\sigma$=' num2str(std_xyz,2)]);
hl = legend('$\pm \sigma$','Ensemble MSD XYZ','Fit points','Fitted model','Location','Best');
set([hx,hy,ht,hl],'interpreter','latex','FontSize',15);

% ensemble MSD x
subplot(2,2,2);
pxx = [tau_vec(1:end-1), fliplr(tau_vec(1:end-1))];
pyy = [ensemble_msd_plus_std_x(1:end-1), fliplr(ensemble_msd_minus_std_x(1:end-1))];
fill(pxx, pyy, 1, 'FaceColor', [0.5,0.5,0.5],'EdgeColor', 'none', 'FaceAlpha',0.2);
hold on;
plot(tau_vec, ensemble_msd_x,'b','Linewidth',2);
hold on;
plot(tau_fit, ensemble_fit_x,'or','Linewidth',2);
hold on;
plot(tau_msd_0, ensemble_msd_0_x,'--k','Linewidth',2);
hold off;  
grid on;
xlim([0, 30]);
hx = xlabel('$\tau \ \left[Sec\right]$');
hy = ylabel('MSD $\left[\frac{\mu m^2}{sec}\right]$');
ht = title(['Ensemble MSD X: D=' num2str(D_x,2) ', $\sigma$=' num2str(std_x,2)]);
hl = legend('$\pm \sigma$','Ensemble MSD X','Fit points','Fitted model','Location','Best');
set([hx,hy,ht,hl],'interpreter','latex','FontSize',15);

% ensemble MSD y
subplot(2,2,3);
pxx = [tau_vec(1:end-1), fliplr(tau_vec(1:end-1))];
pyy = [ensemble_msd_plus_std_y(1:end-1), fliplr(ensemble_msd_minus_std_y(1:end-1))];
fill(pxx, pyy, 1, 'FaceColor', [0.5,0.5,0.5],'EdgeColor', 'none', 'FaceAlpha',0.2);
hold on;
plot(tau_vec, ensemble_msd_y,'b','Linewidth',2);
hold on;
plot(tau_fit, ensemble_fit_y,'or','Linewidth',2);
hold on;
plot(tau_msd_0, ensemble_msd_0_y,'--k','Linewidth',2);
hold off;  
grid on;
xlim([0, 30]);
hx = xlabel('$\tau \ \left[Sec\right]$');
hy = ylabel('MSD $\left[\frac{\mu m^2}{sec}\right]$');
ht = title(['Ensemble MSD Y: D=' num2str(D_y,2) ', $\sigma$=' num2str(std_y,2)]);
hl = legend('$\pm \sigma$','Ensemble MSD Y','Fit points','Fitted model','Location','Best');
set([hx,hy,ht,hl],'interpreter','latex','FontSize',15);

% ensemble MSD z
subplot(2,2,4);
pxx = [tau_vec(1:end-1), fliplr(tau_vec(1:end-1))];
pyy = [ensemble_msd_plus_std_z(1:end-1), fliplr(ensemble_msd_minus_std_z(1:end-1))];
fill(pxx, pyy, 1, 'FaceColor', [0.5,0.5,0.5],'EdgeColor', 'none', 'FaceAlpha',0.2);
hold on;
plot(tau_vec, ensemble_msd_z,'b','Linewidth',2);
hold on;
plot(tau_fit, ensemble_fit_z,'or','Linewidth',2);
hold on;
plot(tau_msd_0, ensemble_msd_0_z,'--k','Linewidth',2);
hold off;  
grid on;
xlim([0, 30]);
hx = xlabel('$\tau \ \left[Sec\right]$');
hy = ylabel('MSD $\left[\frac{\mu m^2}{sec}\right]$');
ht = title(['Ensemble MSD Z: D=' num2str(D_z,2) ', $\sigma$=' num2str(std_z,2)]);
hl = legend('$\pm \sigma$','Ensemble MSD Z','Fit points','Fitted model','Location','Best');
set([hx,hy,ht,hl],'interpreter','latex','FontSize',15);
linkaxes;

% saving everything
save(strcat('ensemble_msd_fits_TP_',str_density,'.mat'),'tau_fit','tau_msd_0',...
    'ensemble_fit_xyz','ensemble_msd_0_xyz','pxyz','std_xyz','D_xyz',...
    'ensemble_fit_x','ensemble_msd_0_x','px','std_x','D_x',...
    'ensemble_fit_y','ensemble_msd_0_y','py','std_y','D_y',...
    'ensemble_fit_z','ensemble_msd_0_z','pz','std_z','D_z');
