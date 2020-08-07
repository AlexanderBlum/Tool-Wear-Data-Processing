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
    [serial, theta, brand] = build_filenames_using_reg_exp_r8(fname);
    tool = strjoin({brand, serial, theta, 'deg'});
else
    tool = 'no struct';
end

end