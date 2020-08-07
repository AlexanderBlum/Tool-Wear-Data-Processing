function data_processing_main_r10(input_data)

save_data = 0;
plot_data = 0;
save_dir = ['C:\Users\Alexa\Google Drive\01 School\02 Research\01 Tool Wear Research\02 Experiment data\',...
            'Tool Wear Experiment Data\01 Good Data\KY 314561 0 deg 20% O2 - Repeat'];
res_axis = 'small';


%% close any uiwaitbar's that were left open
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);

%% input parameters
% GENERAL PARAMETERS     
% parameter list:
% fixed paramaters
text_size = [9.0 5.0];
pixels     = 1024;
fov        = 420;
map_crop   = 100;
dx_interp = .001;

n_plunges = length(input_data);

roughPlunges(n_plunges,1) = SurfAnalysis();
alignedPlunges(n_plunges,1) = SurfAnalysis();
plungeResiduals(n_plunges,1) = SurfAnalysis();
plungeMaps(n_plunges,1) = SurfAnalysis();

%% create file name parts
try
    tool = create_filename(input_data);
catch
    tool = 'unknown (regexp failed)';
end

%% main program
%% prepare data for final alignment
% change z scale from meters to um
for ii = 1:n_plunges
    plungeMaps(ii) = SurfAnalysis(input_data(ii).phaseMap, fov/pixels, 'm');    
    plungeMaps(ii).PhaseMap = plungeMaps(ii).PhaseMap.*10^6;
    plungeMaps(ii).Xscale = 'um';    
    plungeMaps(ii).Zscale = 'um';
end

clear input_data 

% check orientation: take std dev of col avg data along dim 2, if it is greater than 1 the data is aligned
if std(nanmean(plungeMaps(1).PhaseMap, 2)) < 1
    for ii = 1:n_plunges
        plungeMaps(ii).PhaseMap = plungeMaps(ii).PhaseMap';
    end
end

% remove best fit plane
pre_fig = gcf;
pre_ax = gca;     
hold off
set(gcf, 'Units', 'Normalized', 'OuterPosition', [.01 .05 .95 .9]);
for ii = 1:n_plunges
    plot(pre_ax, nanmean(plungeMaps(ii).PhaseMap,2));
    title(pre_ax, 'Remove Plane: Select Flats');
    [ind, ~] = ginput(2);
    ind = floor(ind);
%     close(gcf);   
    ref_phasemap = plungeMaps(ii).PhaseMap;  
    ref_phasemap(ind(1):ind(2),:) = nan;
    plane = fit_plane( ref_phasemap );
    plungeMaps(ii).PhaseMap = plungeMaps(ii).PhaseMap - plane;
end

clear ref_phasemap plane

% rotate
    axes(pre_ax);
for ii = 1:n_plunges
    Plot(plungeMaps(ii));
    title(pre_ax, 'Choose slice locations from left to right');
    [ind, ~] = ginput(2); %ind(1) is front, ind(2) is back
    ind = floor(ind);
    
    fs = plungeMaps(ii).PhaseMap(:,ind(1));
    bs = plungeMaps(ii).PhaseMap(:,ind(2));
    
    pixel_distance = diff(ind);  

    plot(pre_ax, fs);
    title(pre_ax, 'Rotation Front Slice');

    [ind, ~] = ginput(1);
    front_ind = floor(ind);

    plot(pre_ax, bs);
    title('Rotation Back Slice');

    [ind, ~] = ginput(1);
    back_ind = floor(ind);

    rotation_angle = atand((front_ind-back_ind)/pixel_distance);    
    plungeMaps(ii).PhaseMap = imrotate(plungeMaps(ii).PhaseMap, -rotation_angle, 'bilinear', 'crop'); 

    plungeMaps(ii).PhaseMap = plungeMaps(ii).PhaseMap(map_crop:end-map_crop,map_crop:end-map_crop);
end
close(pre_fig);
clear pixel_distance rotation_angle front_ind back_ind ind

% column average and interp to finer spacing
for ii = 1:n_plunges
    roughPlunges(ii).dx = fov/pixels;
    roughPlunges(ii).Xscale = 'um';
    roughPlunges(ii).Zscale = 'um';    
    roughPlunges(ii).Trace = nanmean(plungeMaps(ii).PhaseMap, 2);
    roughPlunges(ii).Trace = interp1(roughPlunges(ii).X                                                                    ,... % old x
                                roughPlunges(ii).Trace                                                            ,... % z
                                linspace(0, floor(roughPlunges(ii).X(end)/dx_interp)*dx_interp, floor(roughPlunges(ii).X(end)/dx_interp))',... % x interp
                                'linear');
    roughPlunges(ii).dx = dx_interp;
end

clear dx_interp

%% perform final alignment
alignedPlunges = align_plunges_trim(roughPlunges);
for ii = 1:n_plunges
    alignedPlunges(ii).dx = roughPlunges(ii).dx;
    alignedPlunges(ii).Xscale = 'um';
    alignedPlunges(ii).Zscale = 'um';        
end

% create fit residual for zero plunge
% this function needs to received a TRIMMED plunge
z0_fit = fit_plunge(alignedPlunges(1).dx, alignedPlunges(1).Trace, 2);

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
    plungeResiduals(ii).Trace = alignedPlunges(ii).Trace - z0_fit;
    %% STOPPED HERE LAST NIGHT. NEED TO CHANGE RESIDUAL FROM CELL ARRAY TO OBJECT IN ALL CODE BELOW THIS
end

%% plots

x = linspace(0, p_length*dx_interp, p_length);
x = x - mean(x);

clear fs bs ii ind map_crop n_plunges p_length map pixels fov
%% plot data
if plot_data == 1
    figs = plot_results(alignedPlunges, residual, res_axis, text_size, 0, 0);
end
%% save data
if save_data == 1 
    % save processed data to .mat file
    % need to fix this so i am properly saving processed x and z data
    save([save_dir, tool, ' processed SWLI data.mat'],...
         'maps', 'plunges', 'alignedPlunges', 'residual');
end
