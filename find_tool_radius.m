function [r, fitresult, debug] = find_tool_radius(x, z)
%CREATEFIT(TEST_PLUNGE_X,TEST_PLUNGE)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : test_plunge_x
%      Y Output: test_plunge
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( x, z );

% Set up fittype and options.
ft = fittype( '-sqrt(r^2-(x-h)^2)+k', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-150 400 400];
opts.StartPoint = [215 500 500];
opts.Upper = [300 600 600];

% Fit model to data.
[fitresult, ~] = fit( xData, yData, ft, opts );

end