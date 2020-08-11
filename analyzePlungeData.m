
zScaleFactor = 10^6;
userSelect = 0;
analyzeData = PlungeProcessing(output(1:3), 0);

%% main program
%% prepare data for final alignment
% change z scale from meters to um
for ii = 1:analyzeData.Nplunges
    analyzeData.PlungeMapData(ii).PhaseMap = analyzeData.PlungeMapData(ii).PhaseMap*zScaleFactor;
    analyzeData.PlungeMapData(ii).Zscale   = 'um';
end

analyzeData.CheckMapOrientation();

analyzeData.RemovePlane();

analyzeData.RotateMaps();

for ii = 1:analyzeData.Nplunges
    analyzeData.PlungeTraceData(ii).dx = analyzeData.PlungeMapData(ii).dx;
    analyzeData.PlungeTraceData(ii).Xscale = analyzeData.PlungeMapData(ii).Xscale;
    analyzeData.PlungeTraceData(ii).Zscale = analyzeData.PlungeMapData(ii).Zscale;
end

analyzeData.GetAvgSlices();

for ii = 1:analyzeData.Nplunges
    analyzeData.PlungeTraceData(ii).InterpTrace(obj.dxInterp);
end % interp plunge trace data

analyzeData.InterpPlungeTraces();

analyzeData.AlignPlungesTrim();

analyzeData.FitFirstPlunge();

analyzeData.CalcResiduals();

analyzeData.PlotResults();
