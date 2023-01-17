function [SP1, SP2] = produceNumberOfSurfaceProtein(xmin, xmax, TotalNumberOfSurfaceProtein)
% Calculating number of HA/NA depending on total surface proteins
% and two ratios observed by Harris et al.

x = xmin+rand(1)*(xmax-xmin); % using uniform distribution across the ratios

% HA present:
SP1 = round((TotalNumberOfSurfaceProtein/(x+1))*x);

% NA present:
SP2 = round((TotalNumberOfSurfaceProtein/(x+1))*1);
