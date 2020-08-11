function [serial, theta, brand] = build_filenames_using_reg_exp(str)

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

end