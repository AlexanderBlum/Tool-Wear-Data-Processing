function plot_debug_data(plunges)

%% create figures
h.unaligned_plunges = figure('Visible'      , 'Off'              ,...
                             'Units'        , 'Normalized'       ,...
                             'NumberTitle'  , 'off'              ,...
                             'Name'         , 'Unaligned plunges',...
                             'OuterPosition', [  0  .5  .3  .45 ]);  
                         
h.col_avg_std  =      figure('Visible'      , 'Off'           ,...
                             'Units'        , 'Normalized'    ,...
                             'NumberTitle'  , 'off'           ,...
                             'Name'         , 'Col avg std dev'  ,...
                             'OuterPosition', [ .6   .05  .3  .45 ]);                             
                             
%% Plot unaligned plunges                         
    unaligned_plunge_ax = axes(h.unaligned_plunges); hold(unaligned_plunge_ax, 'on');
    set(h.unaligned_plunges,'CurrentAxes', unaligned_plunge_ax);
    offset = 0;
    for ii = 1:length(plunges)
        plot(unaligned_plunge_ax, nanmean(plunges(ii).PhaseMap, 2) + offset);
        offset = offset + .5;
    end
%     cellfun(@plot, x_interp, z_interp); 
    xlabel(unaligned_plunge_ax, 'Plunge position (\mum)');
    ylabel(unaligned_plunge_ax, 'Plunge depth (\mum)');

    legend(unaligned_plunge_ax, {'plunge 0', 'plunge 1', 'plunge 2', 'plunge 3', ...
                  'plunge 4', 'plunge 5'}                       , ...
                  'Location', 'best', 'AutoUpdate', 'off');
    
%% plot std deviation of col avg plunges    
    col_avg_std_ax = axes(h.col_avg_std); hold(col_avg_std_ax, 'on');
    set(h.col_avg_std,'CurrentAxes', col_avg_std_ax);
    offset = 0;
    for ii = 1:length(plunges)
        plot(col_avg_std_ax, nanstd(plunges(ii).PhaseMap,1, 2) + offset);
        offset = offset + .25;
    end

    xlabel(col_avg_std_ax, 'Plunge position (index)');
    ylabel(col_avg_std_ax, 'Col avg std dev (\mum)');

%% display plots

hn = fieldnames(h); % create cell array of all field names in structure
for ii = 1:numel(hn)
    % loop structure by field name, set Units and NumberTitle properties
    h.(hn{ii}).Visible = 'On';  
end
end