
zScaleFactor = 10^6;
userSelect = 0;
analyzeData = PlungeProcessing(output, 0);

%% main program
%% prepare data for final alignment
% change z scale from meters to um
for ii = 1:analyzeData.Nplunges
    analyzeData.PlungeData(ii).PhaseMap = analyzeData.PlungeData(ii).PhaseMap*zScaleFactor;
    analyzeData.PlungeData(ii).Zscale   = 'um';
end

analyzeData.CheckMapOrientation();

analyzeData.RemovePlane()

analyzeData.RotateMaps()

analyzeData.