classdef PlungeProcessing
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        PlungeMaps
        PlungeTraces
        Name
    end
    
    methods
        function obj = PlungeProcessing(input_data, n_plunges, dx, map_crop, dx_interp)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            plunges(n_plunges,1) = SurfAnalysis();
            obj.PlungeMaps = plunges;
            obj.PlungeTraces = plunges;
            obj.Name = create_filename(input_data);           
        end
        
        function checkMapOrientation(obj)
        end
        
        function removeBestFitPlane(obj)
        end
        
        function rotate
    end
end

%% helper functions for the class methods
function tool = create_filename(input_data)
% 
% if isstruct(input_data)
%     fname = input_data(1).name(1:end-9) ; % get filename from input structure
%     % parse fname for the components we care about
%     [serial, theta, brand] = build_filenames_using_reg_exp(fname);
%     tool = strjoin({brand, serial, theta, 'deg'});
% else
%     tool = 'no struct';
% end

% updated vesion , July 2019. splits the string using a delimiter at 'deg'
if isstruct(input_data)
    fname = strsplit(input_data(1).name, ' deg') ; % get filename from input structure
    fname = [fname{1} ' deg'];
    % parse fname for the components we care about
    [serial, theta, brand] = build_filenames_using_reg_exp(fname);
    tool = strjoin({brand, serial, theta, 'deg'});
else
    tool = 'no struct';
end

end

