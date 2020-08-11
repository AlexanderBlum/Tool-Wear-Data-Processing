classdef PlungeProcessing < handle
    
    properties
        PlungeMapData   SurfAnalysis
        PlungeTraceData SurfAnalysis
        PlungeResidData SurfAnalysis
        Name
        Nplunges
        
    end
    properties (Constant)
        dxInterp  = 0.001; % micrometers
        MapPixels = 1024;
        MapFOV    = 420; % micrometers
        MapCrop   = 100;
        R         = 500; % tool radius in micrometers
        % deg, from maxSlope = atand(sqrt(2*doc/(R-2*doc)));
        SlopeTol  = 10;      
    end
    properties (Access = private)
        UserSelect = 0;
        EdgeFindStartDist = 20; % in micrometers
        EdgeFindStartSide {mustBeMember(EdgeFindStartSide,{'l','r'})} = 'l';
    end
    
    methods
        function obj = PlungeProcessing(inputData, userSelect)
        % constructor function
            obj.UserSelect = userSelect;
            obj.Nplunges = length(inputData);
            for ii = 1:obj.Nplunges
                MapFOV = obj.MapFOV;
                MapPixels = obj.MapPixels;
                obj.PlungeMapData(ii) = ...
                    SurfAnalysis(inputData(ii).phaseMap, MapFOV/MapPixels, 'm', 'um');
                obj.PlungeTraceData = SurfAnalysis();
                obj.
            end
            obj.Name = createFilename(inputData);   
        end % constructor
 
        function obj = CheckMapOrientation(obj)
            for ii = 1:obj.Nplunges
                if std(nanmean(obj.PlungeMapData(ii).PhaseMap, 2)) < 1
                    obj.PlungeMapData(ii).PhaseMap = obj.PlungeMapData(ii).PhaseMap';
                end
            end
        end % check map orientation
        
        function RemovePlane(obj)
            for ii = 1:obj.Nplunges
                avgSlice = obj.PlungeMapData(ii).GetAvgSlice('col');
                % remove slope and dc bias
                % chunk of slice, the array values are in micrometers
                ind = floor([10, 50]./obj.PlungeMapData(ii).dx); 
                avgSlice = removeTilt(avgSlice, ind);                
                if obj.UserSelect == 0                                    
                    obj.EdgeFindStartSide = 'l';
                    efStartInd = floor(obj.EdgeFindStartDist./obj.PlungeMapData(ii).dx);
                    leftEdgeInd = edgeFinder(avgSlice,...
                                     obj.EdgeFindStartSide,...
                                     efStartInd,...
                                     obj.SlopeTol,...
                                     obj.PlungeMapData(ii).dx);         
                elseif obj.UserSelect == 1
                    if ii == 1
                        figure('Units', 'Normalized', 'OuterPosition', [.01 .05 .95 .9]);
                    end
                    plot(avgSlice)
                    title('Remove Plane: Select Flats');
                    leftEdgeInd = floor(ginput(1));
                    leftEdgeInd = leftEdgeInd(1);
                    if ii == obj.Nplunges
                        close(gcf);  
                    end                    
                end
            % instead of finding both edges, which is buggy
            % find the plunge width with sag equation
            % then add that to the left index to find the right index
                plungeWidth = ...
                    floor(sqrt(abs(min(avgSlice))*8*obj.R)/obj.PlungeMapData(ii).dx);
                indOffset = floor(.05*length(avgSlice));                   
                maskInd = floor([leftEdgeInd - indOffset,...
                                leftEdgeInd + plungeWidth + indOffset]);
                refSurf = SurfAnalysis();
                refSurf.PhaseMap = obj.PlungeMapData(ii).PhaseMap;   
                refSurf.PhaseMap(maskInd(1):maskInd(2), :) = nan;   
                fitPlane = refSurf.FitPlane();
                obj.PlungeMapData(ii).PhaseMap = ...
                    obj.PlungeMapData(ii).PhaseMap - fitPlane;                  
            end
        end % remove best fit plane
    
        function obj = InterpPlungeMaps(obj)
            for ii = 1:obj.Nplunges
                obj.PlungeMapData(ii).InterpMap(obj.dxInterp);
            end
        end % interp plunge maps
        
        function obj = InterpPlungeTraces(obj)
            for ii = 1:obj.Nplunges
                obj.PlungeMapData(ii).InterpTrace(obj.dxInterp);
            end
        end     
        
        function obj = RotateMaps(obj)
        end
        
        function obj = GetAvgSlices(obj)
            for ii = 1:obj.Nplunges
                obj.PlungeMapData(ii).Trace = ...
                    obj.PlungeMapData(ii).GetAvgSlice('col');
            end
        end
        
        function obj = alignPlungesTrim(obj)
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
    start = startIndex              ;
elseif startSide == 'r'
    start = length(z) - start_index  ;
    step  = -step                    ;
end

dz = nan(1, length(z));
for ii = 3:length(z)-2
    % 2nd order accurate finite difference for first derivative
    dz(ii) = atand((-z(ii+2) + 8*(z(ii+1) - z(ii-1)) +z(ii-2))/12/dx);
end

for ii = start:step:floor(length(z)/2)
    deltaZ = dz(ii);
    % difference in z-position from point i to point i + stepSize; if statement
    % below uses this to decide where the edge is based on tol parameter
    if startSide == 'l' && deltaZ < -tol
        ind = ii;
        break
    elseif startSide == 'r' && deltaZ > tol
        ind = ii;
        break     
    end
end   
end

function z = removeTilt(z, ind)
    fs   = z(ind(1):ind(2))                              ; % fit slice
    xFit = linspace(0,length(fs),length(fs))'          ; % x array for fit
    fc   = polyfit(xFit,fs,1)                            ; % fit coefficient
    zFit = polyval(fc,linspace(0,length(z),length(z))) ; % fit polynomial
    z    = z - zFit'                                      ; % z with tilt removed
end


% function maxSlope = maxPlungeSlope(R, doc)
%                 MaxSlope = atand(sqrt(2*doc/(R-2*doc)));
% end

