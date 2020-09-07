function [fitresult, gof, zHat] = fitCircleToTool(x, z)

[x, z] = prepareCurveData( x, z );

ft = fittype( '-sqrt(r^2-(x-h)^2)+k', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-150 400 400];
opts.StartPoint = [215 500 500];
opts.Upper = [300 600 600];

[fitresult, gof] = fit( x, z, ft, opts );
zHat = fitresult(x);
end