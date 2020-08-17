classdef PlungeProcessing < handle
    
    properties
        MapData       SurfAnalysis
        TraceData     SurfAnalysis
        ResidualData  SurfAnalysis
        DebugData     ProcessingDebug
        PlungeZeroFit double
        FileName      char
        ResidCalcType {mustBeMember(ResidCalcType,{'subtractPlungeZero','subtractFit'})} = 'subtractFit'
        TrimmedPlungeLengths double
        UserSelect(1,1) int8 = 0;
        TrimZo(1,1) double
        TrimH(1,1) double
        PlungeFitOrder(1,1) double = 2
        thetaShift(1,1) double
    end
    properties (Hidden)
        
    end
    properties (Constant)
        dxInterp  double = 0.1; % micrometers
        Pixels double = 1024;
        FOV    double = 420; % micrometers
        MapCrop   double = 100;
        R         double = 500; % tool radius in micrometers
        % deg, from maxSlope = atand(sqrt(2*doc/(R-2*doc)));
        SlopeTol  double = 10;      
    end
    properties (Dependent)
        PlungeWidths  double
    end
    properties (SetAccess = immutable)
        Nplunges   int8   
    end
    properties (Access = private)
        RotCrop           = 100;
        EdgeFindStartDist = 20; % in micrometers
        EdgeFindStartSide {mustBeMember(EdgeFindStartSide,{'l','r'})} = 'l';
    end
    
    methods
        function obj = PlungeProcessing(inputData)
        % constructor function
        obj.Nplunges = length(inputData);        
        if isstruct(inputData)
            for ii = 1:obj.Nplunges
                FOV = obj.FOV;
                Pixels = obj.Pixels;
                obj.MapData(ii) = ...
                    SurfAnalysis(inputData(ii).phaseMap, FOV/Pixels, 'm', 'um');
                obj.TraceData = SurfAnalysis();
                obj.ResidualData = SurfAnalysis();
            end
            obj.FileName = createFilename(inputData);              
        else
            obj.MapData = inputData;
            obj.FileName = inputData(1).Name;
            for ii = 1:obj.Nplunges
                obj.TraceData = SurfAnalysis();                
                obj.ResidualData = SurfAnalysis();                
            end
        end
            obj.TrimmedPlungeLengths = nan(1, length(inputData));
            obj.DebugData = ProcessingDebug(length(inputData));
        end % constructor
 
        function obj = CheckMapOrientation(obj)
            for ii = 1:obj.Nplunges
                if std(nanmean(obj.MapData(ii).PhaseMap, 2)) < 1
                    obj.MapData(ii).PhaseMap = obj.MapData(ii).PhaseMap';
                end
            end
        end % check all map orientations
        
        function obj = RemovePlane(obj)
            for ii = 1:obj.Nplunges
                avgSlice = obj.MapData(ii).GetAvgSlice('col');
                % remove slope and dc bias
                % chunk of slice, the array values are in micrometers
                ind = round([10, 50]./obj.MapData(ii).dx); 
                avgSlice = removeTilt(avgSlice, ind);                
                if obj.UserSelect == 0                                    
                    obj.EdgeFindStartSide = 'l';
                    efStartInd = round(obj.EdgeFindStartDist./obj.MapData(ii).dx);
                    leftEdgeInd = edgeFinder(avgSlice,...
                                     obj.EdgeFindStartSide,...
                                     efStartInd,...
                                     obj.SlopeTol,...
                                     obj.MapData(ii).dx);         
                elseif obj.UserSelect == 1
                    if ii == 1
                        figure('Units', 'Normalized', 'OuterPosition', [.01 .05 .95 .9]);
                    end
                    plot(avgSlice)
                    title('Remove Plane: Select Flats');
                    leftEdgeInd = round(ginput(1));
                    leftEdgeInd = leftEdgeInd(1);          
                end
                % instead of finding both edges, which is buggy
                % find the plunge width with sag equation
                % then add that to the left index to find the right index
                plungeWidth = ...
                    round(sqrt(abs(min(avgSlice))*8*obj.R)/obj.MapData(ii).dx);
                indOffset = round(.05*length(avgSlice));                   
                maskInd = round([leftEdgeInd - indOffset,...
                                leftEdgeInd + plungeWidth + indOffset]);
                refSurf = SurfAnalysis();
                refSurf.PhaseMap = obj.MapData(ii).PhaseMap;   
                refSurf.PhaseMap(maskInd(1):maskInd(2), :) = nan;   
                fitPlane = refSurf.FitPlane();
                obj.MapData(ii).PhaseMap = ...
                    obj.MapData(ii).PhaseMap - fitPlane;                  
            end
        end % RemovePlane
    
        function obj = InterpMaps(obj)
            for ii = 1:obj.Nplunges
                obj.MapData(ii).InterpMap(obj.dxInterp);
            end
        end % interp all plunge maps
        
        function obj = InterpTraces(obj)
            for ii = 1:obj.Nplunges
                obj.TraceData(ii).InterpTrace(obj.dxInterp);
            end
        end % interp all plunge traces
        
        function obj = RotateMaps(obj)
            for ii = 1:obj.Nplunges
                % interp         
%                 obj.MapData(1).RotateSurf(-5)                
                frontSliceIndex = round(0.1*obj.MapData(ii).Ncols);
                backSliceIndex  = round(0.9*obj.MapData(ii).Ncols);     
                frontSlice = obj.MapData(ii).GetSlice(frontSliceIndex, 'col');
                backSlice = obj.MapData(ii).GetSlice(backSliceIndex, 'col');      
                slicePixelDist = backSliceIndex - frontSliceIndex;                
                if obj.UserSelect == 0
                    obj.EdgeFindStartSide = 'l';
                    efStartInd = round(obj.EdgeFindStartDist./obj.MapData(ii).dx);                    
                    frontIndex = edgeFinder(frontSlice,...
                                     obj.EdgeFindStartSide,...
                                     efStartInd,...
                                     obj.SlopeTol,...
                                     obj.MapData(ii).dx);  
                    backIndex = edgeFinder(backSlice,...
                                     obj.EdgeFindStartSide,...
                                     efStartInd,...
                                     obj.SlopeTol,...
                                     obj.MapData(ii).dx);                                     
                elseif obj.UserSelect == 1
                    plot(frontSlice)
                    title('Rotation Front Slice');                    
                    frontIndex = round(ginput(1));
                    frontIndex = frontIndex(1);                        

                    plot(backSlice)
                    title('Rotation Back Slice');                    
                    backIndex = round(ginput(1));
                    backIndex = backIndex(1);
                    if ii == obj.Nplunges
                        close(gcf)
                    end
                end
                rotAngle = -atand((frontIndex-backIndex)/slicePixelDist);    
                obj.MapData(ii).RotateSurf(rotAngle);
                obj.MapData(ii).PhaseMap = ...
                    obj.MapData(ii).PhaseMap(obj.RotCrop:end-obj.RotCrop,...
                                                   obj.RotCrop:end-obj.RotCrop);
            end
        end % rotate all maps to common coord
        
        function obj = GetAvgSlices(obj)
            for ii = 1:obj.Nplunges
                obj.TraceData(ii).Trace = ...
                    obj.MapData(ii).GetAvgSlice('col');
            end
        end
        
        function obj = AlignPlungesTrim(obj)
            dzStep = 0.0025;
%             dzStep = 0.0005;
            if obj.UserSelect == 1
                obj.TraceData(1).Plot()
                title('1st click = z0, |2nd click - 3rd click| = h');
                [~, zGinput] = ginput(3);
                close(gcf)
                obj.TrimZo = zGinput(1);     
                obj.TrimH = abs(diff(zGinput(2:3)));                  
            end          
            trimWidth = sqrt(8*obj.TrimH*obj.R);
            for ii = 1:obj.Nplunges
                % trim to common plunge width                
                obj.TraceData(ii).Trace = ...
                    trimPlunge(obj.TraceData(ii).Trace,...
                               trimWidth, dzStep, obj.TraceData(ii).X, obj.TrimZo);
                % remove DC offset
                obj.TraceData(ii).Trace = obj.TraceData(ii).Trace ...
                                        - mean([obj.TraceData(ii).Trace(1),...
                                               obj.TraceData(ii).Trace(end)]);
                obj.TrimmedPlungeLengths(ii) = ...
                    length(obj.TraceData(ii).Trace);
            end            
        end

        function obj = AlignPlungesSlope(obj)                      
        unworn
        
        end
        
        function obj = getThetaShift(obj)
            efStartInd = round(obj.EdgeFindStartDist/obj.TraceData(1).dx);           
            % trim plunges
            leftInd = edgeFinder(obj.TraceData(1).Trace,...
                'l', efStartInd, obj.SlopeTol, obj.TraceData(1).dx);
            rightInd = edgeFinder(obj.TraceData(1).Trace,...
                'r', efStartInd, obj.SlopeTol, obj.TraceData(1).dx);
            z = obj.TraceData(1).Trace(leftInd:rightInd);
            x = obj.TraceData(1).X(leftInd:rightInd);
            x = x - mean(x);
            shiftLim = 16;
            shiftStep = 4;            
            obj.thetaShift = calcIndexShiftVsTheta(x, z, shiftLim, shiftStep);                
        end
        
        function obj = FitFirstPlunge(obj, n)
            PlungeFit = polyfit(obj.TraceData(1).X',...
                                obj.TraceData(1).Trace,...
                                n);
            obj.PlungeZeroFit = polyval(PlungeFit, obj.TraceData(1).X)';
        end       
  
        function obj = CalcResiduals(obj)
            % find shortest z vector. truncate all z and x vectors to that length.
            % remove DC offset from z vectors. shift x vectors to be symmetric
            % should I reinterp all of the vectors to same length instead?
%             maxDelta = max(obj.TrimmedPlungeLengths) - min(obj.TrimmedPlungeLengths);
            minLength = min(obj.TrimmedPlungeLengths);
                switch obj.ResidCalcType
                    case 'subtractFit'
                        for ii = 1:obj.Nplunges
                            obj.ResidualData(ii).Trace = obj.TraceData(ii).Trace(1:minLength) ...
                                                       - obj.PlungeZeroFit(1:minLength);
                        end
                    case 'subtractPlungeZero'
                        for ii = 2:obj.Nplunges
                            obj.ResidualData(ii).Trace = obj.TraceData(ii).Trace(1:minLength) ...
                                                       - obj.TraceData(1).Trace(1:minLength);
                        end
                end            
        end
        
        function obj = PlotResults(obj)
        end
        
        function obj = PlotSlices(obj)
            figure;
            hold on;
            for ii = 1:obj.Nplunges
                obj.TraceData(ii).Plot();
            end
%             axis([-
        end
        
        function obj = PlotResiduals(obj)
            figure;
            hold on;
            for ii = 1:obj.Nplunges
                obj.ResidualData(ii).Plot();
            end      
        end
        
        function PlungeWidths = get.PlungeWidths(obj)
            PlungeWidths = nan(1, obj.Nplunges);
            for ii = 1:obj.Nplunges                
                obj.EdgeFindStartSide = 'l';
                efStartInd = round(obj.EdgeFindStartDist./obj.TraceData(ii).dx);
                leftEdgeIndex = edgeFinder(obj.TraceData(ii).Trace,...
                                           obj.EdgeFindStartSide,...
                                           efStartInd,...
                                           obj.SlopeTol,...
                                           obj.TraceData(ii).dx);
                obj.EdgeFindStartSide = 'r';                             
                rightEdgeIndex = edgeFinder(obj.TraceData(ii).Trace,...
                                           obj.EdgeFindStartSide,...
                                           efStartInd,...
                                           obj.SlopeTol,...
                                           obj.TraceData(ii).dx);
                PlungeWidths(ii) = abs(leftEdgeIndex-rightEdgeIndex);
            end
        end
        
        function obj = ProcessPlunges(obj)
            for ii = 1:obj.Nplunges
                obj.MapData(ii).PhaseMap = obj.MapData(ii).PhaseMap;
                obj.MapData(ii).Zscale   = 'um';
            end

            obj.CheckMapOrientation();
            obj.InterpMaps()
%             obj.RemovePlane();

            obj.RotateMaps();

            obj.GetAvgSlices();

            % add most current property data to the plunge traces
            for ii = 1:obj.Nplunges
                obj.TraceData(ii).dx = obj.MapData(ii).dx;
                obj.TraceData(ii).Xscale = obj.MapData(ii).Xscale;
                obj.TraceData(ii).Zscale = obj.MapData(ii).Zscale;
            end

%             obj.InterpTraces();

            % obj.PlotSlices();   % see if column average went well

            obj.AlignPlungesTrim();
%             obj = getThetaShift(obj);
            
%             obj.AlignPlungesSlope();

%             obj.PlotSlices();   % see if alignment went well

%             obj.FitFirstPlunge(obj.PlungeFitOrder);

            % add most current property data to the plunge residuals
            for ii = 1:obj.Nplunges
                obj.ResidualData(ii).dx = obj.MapData(ii).dx;
                obj.ResidualData(ii).Xscale = obj.MapData(ii).Xscale;
                obj.ResidualData(ii).Zscale = obj.MapData(ii).Zscale;
            end

        obj.CalcResiduals();            
        end
    end

end

%% helper functions for the class methods
function fname = createFilename(input_data)

% updated vesion , July 2019. splits the string using a delimiter at 'deg'
if isstruct(input_data)
    fname = strsplit(input_data(1).name, ' deg') ; % get filename from input structure
    fname = [fname{1} ' deg'];
    % parse fname for the components we care about
else
    fname = 'no struct';
end

[serial, theta, brand] = findFilenameParts(fname);
fname = strjoin({brand, serial, theta, 'deg'});
    
end % createFilename

function [serial, theta, brand] = findFilenameParts(str)

% find serial number
if regexp(str, '[0-9]{6,7}|unknown')
    serial = regexp(str, '[0-9]{6,7}|unknown', 'match');
else
    serial = 'serial not found';
end

% find angle experiment was performed at
if regexp(str, '\s[-]?[0-9]{1,3}([.]|pt)?[0-9]?\s')
    theta = regexp(str, '\s[-]?[0-9]{1,3}([.]|pt)?[0-9]?\s', 'match');
    theta = strrep(theta{1}(2:end-1), 'pt', '.');
else
    theta = 'unknown';
end

% find brand of tool
if regexp(str, 'K.*Y')
    brand = 'KY';
elseif regexp(str, 'ET|E.*h')
    brand = 'ET';
else
    brand = 'unknown';
end

serial = char(serial);
brand = char(brand);
theta = char(theta);

if regexp(str, 'fake data')
    serial = 'fake';
    brand  = 'fake';
    theta  = 'fake';
end

end % find filename parts

function [ ind ] = edgeFinder(z, startSide, startIndex, tol, dx)
% This function finds the left or right hand edge of any given
% "slice" of the PhaseMap data by looking for a sharp dropoff
% in z-position

step = 1;

% logic for starting the loop, based on side of plunge we want
if  startSide == 'l'
    start = startIndex            ;
elseif startSide == 'r'
    start = length(z) - startIndex;
    step  = -step                 ;
end

dz = nan(1, length(z));
for ii = 3:length(z)-2
    % 2nd order accurate finite difference for first derivative
    dz(ii) = atand((-z(ii+2) + 8*(z(ii+1) - z(ii-1)) +z(ii-2))/12/dx);
end

for ii = start:step:round(length(z)/2)
    deltaZ = dz(ii);
    % difference in z-position from point i to point i + stepSize; if statement
    % below uses this to decide where the edge is based on tol parameter
    if startSide == 'l' && deltaZ < -tol
        ind = ii;
        return
    elseif startSide == 'r' && deltaZ > tol
        ind = ii;
        return     
    end
end
error('No edge found; max deltaZ found is %.2f, min is %.2f', max(dz), min(dz));
end

function z = removeTilt(z, ind)
    fs   = z(ind(1):ind(2))                              ; % fit slice
    xFit = linspace(0,length(fs),length(fs))'          ; % x array for fit
    fc   = polyfit(xFit,fs,1)                            ; % fit coefficient
    zFit = polyval(fc,linspace(0,length(z),length(z))) ; % fit polynomial
    z    = z - zFit'                                      ; % z with tilt removed
end

function zPlunge = trimPlunge(zPlunge, trimWidth, dzStep, x, z0)

% Function attempts to find f^2/8r location that is the same for all 
% plunges.

% h  - height location to try and locate f at
% dz - amount to step down each loop iteration
% x  - interpolated x data from main script
% z  - interpolated z data from main script
% z0 - position on z axis to begin scanning at

for ii = 1:round(length(x)/2)
    % find first and last values in x array that correspond
    % to where z is below the current scan value
    xPair = [x(find(zPlunge < z0, 1, 'first'))  x(find(zPlunge < z0, 1, 'last'))]; 
    xDist = max(xPair) - min(xPair);
    % if only one value that meets the criteria is found, move down in z
    if length(xPair) == 1
        z0 = z0 - dzStep;
    % round x_dist and f to two decimal places. if x_dist is less than or
    % equal to f, then trim at this location.
    elseif round(xDist, 4) <= round(trimWidth, 4)
        zPlunge = zPlunge(x >= min(xPair) & x <= max(xPair));
        break
    end
    z0 = z0 - dzStep;
end
end

function thetaShift = calcIndexShiftVsTheta(x, z, shiftLim, shiftStep)
z = z';
% calculate relationship between resid
% slope and index shift
shifts = -shiftLim:shiftStep:shiftLim;
shiftResidual = nan(length(x), length(shifts)); 
zShift = nan(length(x), length(shifts)); 
pFit = nan(length(shifts), 2);

leftTrimInd = round(0.2*length(x));
rightTrimInd  = round(0.8*length(x)); 

for ii = 1:length(shifts)
    zShift(:,ii) = circshift(z, shifts(ii));
    shiftResidual(:,ii) = zShift(:,ii) - z;
    pFit(ii, :) = polyfit(x(leftTrimInd:rightTrimInd),...
                          shiftResidual(leftTrimInd:rightTrimInd, ii), 1);
end
thetaShift = polyfit(shifts', pFit(:, 1), 1);
thetaShift = thetaShift(1);
end
% function maxSlope = maxPlungeSlope(R, doc)
%                 MaxSlope = atand(sqrt(2*doc/(R-2*doc)));
% end

