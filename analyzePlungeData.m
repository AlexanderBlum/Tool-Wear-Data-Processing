
zScaleFactor = 10^6;
userSelect = 0;
processedData = PlungeProcessing(output(1:3));
processedData.ResidCalcType = 'subtractPlungeZero';
processedData.UserSelect = 0;
processedData.TrimH = 7;
processedData.TrimZo = -2;
% ResidCalcType = 'subtractPlungeZero';

%% main program
%% prepare data for final alignment
% change z scale from meters to um
for ii = 1:processedData.Nplunges
    processedData.MapData(ii).PhaseMap = processedData.MapData(ii).PhaseMap*zScaleFactor;
    processedData.MapData(ii).Zscale   = 'um';
end

processedData.CheckMapOrientation();

processedData.RemovePlane();

processedData.RotateMaps();

processedData.GetAvgSlices();

% add most current property data to the plunge traces
for ii = 1:processedData.Nplunges
    processedData.TraceData(ii).dx = processedData.MapData(ii).dx;
    processedData.TraceData(ii).Xscale = processedData.MapData(ii).Xscale;
    processedData.TraceData(ii).Zscale = processedData.MapData(ii).Zscale;
end

processedData.InterpTraces();

% analyzeData.PlotSlices();   % see if column average went well

% processedData.PlungeWidths   % only works on untrimmed plunges

% processedData.AlignPlungesTrim();

processedData.AlignPlungesSlope();

% processedData.PlotSlices();   % see if alignment went well

processedData.FitFirstPlunge(3);

% add most current property data to the plunge residuals
for ii = 1:processedData.Nplunges
    processedData.ResidualData(ii).dx = processedData.MapData(ii).dx;
    processedData.ResidualData(ii).Xscale = processedData.MapData(ii).Xscale;
    processedData.ResidualData(ii).Zscale = processedData.MapData(ii).Zscale;
end

processedData.CalcResiduals();

processedData.PlotResiduals();
% processedData.PlotResults();
