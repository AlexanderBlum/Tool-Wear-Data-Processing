function data_processing_main_r11(input_data)
%% close any uiwaitbar's that were left open
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);

%% input parameters
% GENERAL PARAMETERS     
% parameter list:
% fixed paramaters
pixels     = 1024;
fov        = 420;
map_crop   = 100;
dx_interp = .001;

n_plunges = length(input_data);
maps(n_plunges,1) = SurfAnalysis();
plunges(n_plunges,1) = SurfAnalysis();

plungeProcessingObj = PlungeProcessing(input_data, n_plunges, fov/pixels, dx_interp)

%% create object array
for ii = 1:n_plunges
    plungeProcessingObj.PlungeMaps(ii) = SurfAnalysis(input_data(ii).phaseMap, fov/pixels, 'm');
end

%% main program
%% prepare data for final alignment
% change z scale from meters to um
for ii = 1:n_plunges
    plungeProcessingObj.PlungeMaps(ii) = plungeProcessingObj.PlungeMaps(ii).*10^6;
    maps(ii).Zscale = 'um';
end

% check orientation: take std dev of col avg data along dim 2, if it is greater than 1 the data is aligned
if std(nanmean(plungeProcessingObj.PlungeMaps(1).PhaseMap, 2)) < 1
    for ii = 1:n_plunges
        plungeProcessingObj.PlungeMaps(1).PhaseMap = plungeProcessingObj.PlungeMaps(1).PhaseMap';
    end
end

% remove best fit plane
for ii = 1:n_plunges
    figure('Units', 'Normalized', 'OuterPosition', [.01 .05 .95 .9]);
    plot(nanmean(maps(ii).PhaseMap,2));
    title('Remove Plane: Select Flats');
    [x, ~] = ginput(2);
    x = floor(x);
    close(gcf);   
    ref_phasemap = maps(ii).PhaseMap;  
    ref_phasemap(x(1):x(2),:) = nan;
    plane = fit_plane( ref_phasemap );
    maps(ii).PhaseMap = maps(ii).PhaseMap - plane;
end

% rotate
for ii = 1:n_plunges
    plot(maps(ii));
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [.01 .05 .95 .9]);
    title('Choose slice locations from left to right');
    [ind, ~] = ginput(2); %ind(1) is front, ind(2) is back
    ind = floor(ind);
    
    fs = maps(ii).PhaseMap(:,ind(1));
    bs = maps(ii).PhaseMap(:,ind(2));
    
    pixel_distance = diff(ind);
    
    close(gcf);
    
    figure('Units', 'Normalized', 'OuterPosition', [.01 .05 .95 .9]);    
    plot(fs);
    title('Rotation Front Slice');

    [x, ~] = ginput(1);
    x_front = floor(x);
    close(gcf);
    figure('Units', 'Normalized', 'OuterPosition', [.01 .05 .95 .9]);
    plot(bs);
    title('Rotation Back Slice');

    [x, ~] = ginput(1);
    x_back = floor(x);
    close(gcf);    
    rotation_angle = atand((x_front-x_back)/pixel_distance);    
    maps(ii).PhaseMap = imrotate(maps(ii).PhaseMap, -rotation_angle, 'bilinear', 'crop'); 

    maps(ii).PhaseMap = maps(ii).PhaseMap(map_crop:end-map_crop,map_crop:end-map_crop);
end

% column average and interp to finer spacing
for ii = 1:n_plunges
    plunges(ii).Trace = nanmean(maps(ii).PhaseMap, 2);
    x = linspace(0, length(plunges(ii).Trace)*fov/pixels, length(plunges(ii).Trace))';
    Ninterp = floor(x(end)/dx_interp);
    x_interp = linspace(0, Ninterp*dx_interp, Ninterp)';  
    plunges(ii).Trace = interp1(x, plunges(ii).Trace, x_interp, 'linear');
    plunges(ii).dx = dx_interp;
end

%% perform final alignment
alignedPlunges = align_plunges_trim(plunges);

% create fit residual for zero plunge
% this function needs to received a TRIMMED plunge
z0_fit = fit_plunge(dx_interp, alignedPlunges(1).Trace, 2);

%% get residuals
% find shortest z vector. truncate all z and x vectors to that length.
% remove DC offset from z vectors. shift x vectors to be symmetric
% should I reinterp all of the vectors to same length instead?
p_length = nan(n_plunges,1);
for ii = 1:n_plunges
    p_length(ii) = length(alignedPlunges(ii).Trace);
end
p_length_max_delta = (max(p_length) - min(p_length))*dx_interp; %#ok<NASGU>
p_length = min(p_length);

residual = cell(1,n_plunges);
z0_fit = z0_fit(1:p_length);
for ii = 1:n_plunges
    alignedPlunges(ii).Trace = alignedPlunges(ii).Trace(1:p_length);
    residual{ii} = alignedPlunges(ii).Trace - z0_fit;
end

figure
for ii = 1:n_plunges
    x = linspace(0, length(residual{ii})*dx_interp, length(residual{ii}));
    x = x - mean(x);
    plot(x, residual{ii});
    hold on
end
% axis([-100, 100, -.1, .3]);
axis([-100  100  -0.1  1.5]);
xlabel('Location On Tool Edge (\mum)');
ylabel('Edge Recession (nm)');
title(tool);
