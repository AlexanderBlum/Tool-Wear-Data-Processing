function fig = plot_results(x, alignedPlunges, residual, res_axis, text_size, nose_lims, edge_lim)

%% set axis limits

% xlim  will be -100 to 100 for plunge and residual plots
% ylim will be set manually as a parameter for residual plot
% ylim will be set automatically for recession plot
% x = x';
switch res_axis
    case 'small'
    ax_lims = [-100  100  -0.1  0.3];  
    case 'medium'
    ax_lims = [-100  100  -0.1  1.5]; 
    case 'large'
    ax_lims = [-100  100  -0.1  3.5]; 
    case 'extra large'
    ax_lims = [-100  100  -0.1  5];     
end   

wr_ax_lims = [ax_lims(1) ax_lims(2) ax_lims(3) .5*ax_lims(4)];
% plots the plunge and residual cell arrays
ls = {'-', '--', ':', '-.', '-', '--', '-', '--', ':', '-.', '-', '--', '-', '--', ':', '-.', '-', '--'};
ls = ls(1:length(alignedPlunges));

%% create figures
fig.plunge            = figure('Visible'      , 'off'              ,...
                               'Units'        , 'Normalized'       ,...
                               'NumberTitle'  , 'off'              ,...
                               'Name'         , 'Aligned plunges',...
                               'OuterPosition', [  0  .5  .3  .5 ]);
                        
fig.resid             = figure('Visible'      , 'off'           ,...
                               'Units'        , 'Normalized'    ,...
                               'NumberTitle'  , 'off'           ,...
                               'Name'         , 'Plunge residuals'  ,...
                               'OuterPosition', [ .3  .5  .3  .5 ]);    
                      
fig.recession_scatter = figure('Visible'      , 'off'           ,...
                               'Units'        , 'Normalized'    ,...
                               'NumberTitle'  , 'off'           ,...
                               'Name'         , 'Recession vs cutting distance'  ,...
                               'OuterPosition', [ .3   0  .3  .5 ]);   
                         
fig.wear_rate         = figure('Visible'      , 'off'           ,...
                               'Units'        , 'Normalized'    ,...
                               'NumberTitle'  , 'off'           ,...
                               'Name'         , 'Wear rate per pass'  ,...
                               'OuterPosition', [ .6  .5  .3  .5 ]);     
                         
fig.wr_scatter        = figure('Visible'      , 'off'           ,...
                               'Units'        , 'Normalized'    ,...
                               'NumberTitle'  , 'off'           ,...
                               'Name'         , 'Wear rate vs cutting distance'  ,...
                               'OuterPosition', [ .6  .0  .3  .5 ]);   
                         
fig.z0fit             = figure('Visible'      , 'off'           ,...
                               'Units'        , 'Normalized'    ,...
                               'NumberTitle'  , 'off'           ,...
                               'Name'         , 'Residual of best fit to plunge 0'  ,...
                               'OuterPosition', [ .0  .0  .3  .5 ]);                            

%% create axes
ax.plunge     = axes(fig.plunge           , 'FontSize', text_size(1),...
                                            'XLim', [ax_lims(1) ax_lims(2)]); 
ax.recession  = axes(fig.recession_scatter, 'FontSize', text_size(1)); 
ax.wear_rate  = axes(fig.wear_rate        , 'FontSize', text_size(1)); 
ax.wr_scatter = axes(fig.wr_scatter       , 'FontSize', text_size(1)); 
ax.z0fit      = axes(fig.z0fit            , 'FontSize', text_size(1)); 

% ax.resid      = axes(fig.resid            , 'FontSize', text_size(1)                          ,...
%                                             'Position', [0.130  0.110+.150  0.775  0.815-.150],...
%                                             'XLim'    , ax_lims(1:2)                          ,...
%                                             'YLim'    , ax_lims(3:4)                          ); 
                        
% ax.chip       = axes(fig.resid            , 'FontSize'     , text_size(1)      ,...
%                                             'Position'     , [.13 .11 .775 .15],...
%                                             'Box'          , 'Off'             ,...
%                                             'YAxisLocation', 'Right'           ,...
%                                             'XLim'         , ax_lims(1:2)      ,...
%                                             'YLim'         , [0  0.3]            ); 

ax.resid      = axes(fig.resid            , 'FontSize', text_size(1)                          ,...
                                            'Position', [0.130  0.180+.145  0.7  0.785-.150],...
                                            'XLim'    , ax_lims(1:2)                          ,...
                                            'YLim'    , ax_lims(3:4)                          ); 
                        
ax.chip       = axes(fig.resid            , 'FontSize'     , text_size(1)      ,...
                                            'Position'     , [.13 .18 .7 .135],...
                                            'Box'          , 'Off'             ,...
                                            'YAxisLocation', 'Right'           ,...
                                            'XLim'         , ax_lims(1:2)      ,...
                                            'YLim'         , [0  2.1]            ); 

ax.resid.XAxis.Visible = 'off';

an = fieldnames(ax);
for ii = 1:numel(an)
    hold(ax.(an{ii}), 'on');
end

%% plot data to figures
% aligned plunges
fig.plunge.Visible = 'on';

for ii = 1:length(alignedPlunges)
    plot(ax.plunge, alignedPlunges(ii).X, alignedPlunges(ii).Trace)
end

% residuals
x = alignedPlunges(1).X;
chip = uncut_chip_r8(5, 2, x); % uncut_chip_r8(doc, feed_rev, x_chip)
plot(ax.chip, x(x > 0), chip(x > 0));

for ii = 1:length(alignedPlunges)
    plot(ax.resid, x, residual{ii});
end

% wear rate
% scatter(ax.wear_rate, cut_dist, r.nose_pv, 'rs'); % plot nose recession vs cut dist
% plot(ax.wear_rate, r.nose_fit, 'r')                 ; % add fit line for nose recession
% scatter(ax.wear_rate, cut_dist, r.edge_pv, 'b+'); % plot edge recession vs cut dist
% plot(ax.wear_rate, r.edge_fit, 'b')                 ; % add fit line for edge recession


plot(ax.z0_fit, x, residual{1}); % plot aligned plunges

%% format the figures
%% format aligned plunges figure (PLUNGE)
set(groot, 'CurrentFigure', fig.plunge);
fig.plunge.CurrentAxes = ax.plunge;

% formatting for aligned plunges plot
xlabel(ax.plunge, 'location along edge (\mum)','FontSize', text_size(1));
ylabel(ax.plunge, 'plunge depth (\mum)','FontSize', text_size(1));
legend(ax.plunge, {'plunge 0', 'plunge 1', 'plunge 2', 'plunge 3' , ...
                   'plunge 4','plunge 5'}                        , ...
                   'Location', 'North', 'FontSize', text_size(2));                           

%% format residual figure (RESID)
set(groot, 'CurrentFigure', fig.resid);
fig.resid.CurrentAxes = ax.resid;

% formatting for residuals plot
yticks(ax.chip, [0 2]);
xlabel(ax.chip, 'location along edge (\mum)', 'FontSize', text_size(1));
chip_string = sprintf('chip \nthickness');
ylabel(ax.chip, [chip_string ' (\mum)'], 'FontSize', 8);
ylabel(ax.resid, 'edge recession (\mum)'     , 'FontSize', text_size(1));
legend(ax.resid, {'plunge 1', 'plunge 2', 'plunge 3', 'plunge 4', ...
                  'plunge 5', 'plunge 6'}                       , ...
                  'Location', 'best', 'FontSize', text_size(2),...
                  'AutoUpdate', 'off');
% add translucent boxed showing where nose and leading edge wear is
% calculated from
% use the variables nose lim and edge lim
% top and bottom of box come from ax_lims index 4, 3, respectively
% so the nose wear box is 
% x = [noselim1 noselim2 noselim1 noselim 2]
% y = [axmin    axmin    axmax    axmax] 
% nose_pgon = polyshape([nose_lims(1) nose_lims(2) nose_lims(2) nose_lims(1)],...
%                       [ax_lims(4)   ax_lims(4)   ax_lims(3)   ax_lims(3)  ]);
% plot(ax.resid, nose_pgon)
% edge_pgon = polyshape([edge_lim   ax_lims(2) ax_lims(2) edge_lim],...
%                       [ax_lims(4) ax_lims(4) ax_lims(3) ax_lims(3)]);
% plot(ax.resid, edge_pgon)

%% format recession vs cut distance figure (RECESSION SCATTER)
% set(groot, 'CurrentFigure', fig.recession_scatter);
% fig.recession_scatter.CurrentAxes = ax.recession;
% 
% % formatting for edge recession scatter plot
% axis(ax.recession, [-.05*max(cut_dist)   1.05*max(cut_dist) ax_lims(3) ax_lims(4)]);
% % set(   ax.recession, 'FontSize', text_size(1));
% xlabel(ax.recession, 'cutting distance (m)' , 'FontSize', text_size(1));
% ylabel(ax.recession, 'edge recession (\mum)', 'FontSize', text_size(1));
% legend(ax.recession, {'nose recession', 'linear fit' , ...
%                       'edge recession', 'linear fit'}, ...
%                       'Location', 'Northwest', 'FontSize', text_size(2));
%                   
% rsqu_nose = ['r^2  nose = ', num2str(round(r.nose_gof.rsquare, 2))];
% rsqu_edge = ['r^2  edge = ', num2str(round(r.edge_gof.rsquare, 2))];
% 
% nose_text = text(.75*max(cut_dist), 1*(0 + ax_lims(4))/2, rsqu_nose);              
% edge_text = text(.75*max(cut_dist), 3*(0 + ax_lims(4))/4, rsqu_edge); 
% 
% set(nose_text, 'fontsize'           , text_size(1), ...
%                'verticalalignment'  , 'bottom'    , ... 
%                'horizontalalignment', 'left'     );
%            
% set(edge_text, 'fontsize'           , text_size(1), ...
%                'verticalalignment'  , 'bottom'    , ... 
%                'horizontalalignment', 'left'     );                  

%% format wear rate residual figure
% set(groot, 'CurrentFigure', fig.wear_rate);
% fig.wear_rate.CurrentAxes = ax.wear_rate;
% 
% % formatting for wear rate residual plot
% axis(ax.wear_rate, wr_ax_lims); 
% % xlim(resid_ax, [-100 100]);
% % set(   ax.wear_rate, 'FontSize', text_size(1));
% 
% cellfun(@plot, x(1:length(ls)-1, 1), dw.wear_rate, ls(1:length(ls)-1)'); % plot aligned plunges
% 
% xlabel(ax.wear_rate, 'location on tool edge (\mum)' , 'FontSize', text_size(1));
% ylabel(ax.wear_rate, 'radial wear rate (nm/m)', 'FontSize', text_size(1));
% legend(ax.wear_rate, {'p1-p0', 'p2-p1', 'p3-p2', 'p4-p3','p5-p4'}, ...
%                      'Location', 'NorthWest', 'FontSize', text_size(2));                 
%% format wear rate vs cut distance scatter figure
% set(groot, 'CurrentFigure', fig.wr_scatter);
% fig.wr_scatter.CurrentAxes = ax.wr_scatter;
% 
% scatter(cut_dist(1:end-1), dw.nose_pv, 'rs'); % plot nose recession vs cut dist
% plot(dw.nose_fit, 'r')                 ; % add fit line for nose recession
% scatter(cut_dist(1:end-1), dw.edge_pv, 'b+'); % plot edge recession vs cut dist
% plot(dw.edge_fit, 'b')                 ; % add fit line for edge recession
% 
% % formatting for edge recession scatter plot
% axis(ax.wr_scatter, [-.05*cut_dist(end-1)   1.05*cut_dist(end-1) ax_lims(3) wr_ax_lims(4)]);
% % set(   ax.wr_scatter, 'FontSize', text_size(1));
% xlabel(ax.wr_scatter, 'cutting distance (m)' , 'FontSize', text_size(1));
% ylabel(ax.wr_scatter, 'wear rate (nm/m)', 'FontSize', text_size(1));
% 
% rsqu_nose = ['r^2  nose = ', num2str(round(dw.nose_gof.rsquare, 2))];
% rsqu_edge = ['r^2  edge = ', num2str(round(dw.edge_gof.rsquare, 2))];
% 
% nose_text = text(.75*max(cut_dist(end-1)), 1*(0 + wr_ax_lims(4))/2, rsqu_nose);              
% edge_text = text(.75*max(cut_dist(end-1)), 3*(0 + wr_ax_lims(4))/4, rsqu_edge); 
% 
% set(nose_text, 'fontsize'           , text_size(1), ...
%                'verticalalignment'  , 'bottom'    , ... 
%                'horizontalalignment', 'left'     );
%            
% set(edge_text, 'fontsize'           , text_size(1), ...
%                'verticalalignment'  , 'bottom'    , ... 
%                'horizontalalignment', 'left'     ); 
%            
% legend(ax.wr_scatter, {'nose wear rate', 'linear fit' , ...
%           'edge wear rate', 'linear fit'}, ...
%           'Location', 'Best', 'FontSize', text_size(2));

%% format figure for zeroth plunge fit
set(groot, 'CurrentFigure', fig.z0fit);
fig.z0fit.CurrentAxes = ax.z0fit;

% formatting for aligned plunges plot
% axis(plunge_ax, ax_lims(1,:)); 
axis(ax.z0fit, 'auto');
xlim(ax.z0fit, [ax_lims(1) ax_lims(2)]);
ylim([-.4  .4]);
% set(ax.z0fit, 'FontSize', text_size(1));
xlabel(ax.z0fit, 'location along edge (\mum)','FontSize', text_size(1));
ylabel(ax.z0fit, 'deviation from best fit (\mum)','FontSize', text_size(1));


%% scale and position the figures for display
fn = fieldnames(fig); % create cell array of all field names in structure
for ii = 1:numel(fn)
    % loop structure by field name, set Units and NumberTitle properties
    fig.(fn{ii}).Units = 'normalized';  
    fig.(fn{ii}).NumberTitle = 'off';
end

for ii = 1:numel(fn)
    % loop structure by field name, set Units and NumberTitle properties
    fig.(fn{ii}).Visible = 'On';  
end
         
end