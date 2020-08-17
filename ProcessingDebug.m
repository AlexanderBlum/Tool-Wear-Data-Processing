classdef ProcessingDebug
    
    properties
        iiCheckOrientations
        
        iiRemovePlanes
        PlaneSlices
        PlaneMaskX
        PlaneMaskZ
        
        iiRotateMaps
        RotateSlices
        RotateLines
        RotateLinesFit
        RotateIndex
        RotateTheta
        
        iiGetAvgSlices
        ColAvgStdDev
        EstPlungeDepth
        
        iiInterpTraces
        iiPlungeAlign
        TrimmedPlungeLengths
        % edge finder
        iiEdgeFinder  
        
        iiCalcResidual
        
    end
    
    methods
        function obj = ProcessingDebug(N)
            obj.PlaneSlices = cell(1,N);
            obj.PlaneMaskX = nan(2, N);
            obj.PlaneMaskZ = nan(2, N);
            obj.RotateSlices = cell(2, N);
            obj.RotateLines = cell(2, N);
            obj.RotateLinesFit = cell(2, N);
            obj.RotateIndex = nan(2, N);
            obj.RotateTheta = nan(1, N);  
            obj.ColAvgStdDev = cell(1, N);     
            obg.EstPlungeDepth = nan(1, N);
        end        
    end
end

