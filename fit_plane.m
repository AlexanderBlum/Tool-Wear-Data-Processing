function [ plane ] = fit_plane( phasemap )
% This function outputs a best fit plane given areal z data 

[x, y] = meshgrid(1:1:size(phasemap,2), 1:1:size(phasemap,1));
X = x(~isnan(phasemap));
Y = y(~isnan(phasemap));
Z = phasemap(~isnan(phasemap));

A = [ones(size(X)) X Y]\Z; % model
plane = ones(size(phasemap)).*A(1) + x.*A(2) + y.*A(3);
end

