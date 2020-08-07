function alignedPlunges = align_plunges_trim(plunges)

% Function attempts to find f^2/8r location that is the same for all 
% plunges.

% h  - height location to try and locate f at
% dz - amount to step down each loop iteration
% x  - interpolated x data from main script
% z  - interpolated z data from main script
% z0 - position on z axis to begin scanning at
dz = 0.0025;

alignedPlunges(length(plunges), 1) = SurfAnalysis();

Plot(plunges(1))
title('1st click = z0, |2nd click - 3rd click| = h');
[~, z_ginput] = ginput(3);
close(gcf)
h = abs(diff(z_ginput(2:3)));

R = 500;
f = sqrt(8*h*R);

for kk = 1:length(plunges)
    z0 = z_ginput(1);
    x = linspace(0, length(plunges(kk).Trace)*plunges(kk).dx, length(plunges(kk).Trace))';
    z = plunges(kk).Trace;
for ii = 1:floor(length(x)/2)
    % find first and last values in x array that correspond
%     to where z is below the current scan value
    x_pair = [x(find(z < z0, 1, 'first'))  x(find(z < z0, 1, 'last'))]; 
    x_dist = max(x_pair) - min(x_pair);
    % if only one value that meets the criteria is found, move down in z
    if length(x_pair) == 1
        z0 = z0 - dz;
    % round x_dist and f to two decimal places. if x_dist is less than or
    % equal to f, then trim at this location.
    elseif round(x_dist, 4) <= round(f, 4)
        z = z(x >= min(x_pair) & x <= max(x_pair));
        alignedPlunges(kk).Trace = z - z(1);
        break
    end
    z0 = z0 - dz;
end
end
end
