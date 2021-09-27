function [x,y,z, n] = obtainProbePositions(labellingAccuracy, numberOfSurfaceProteins, xPoint, yPoint, zPoint);
% Attaching probe onto surface protein for STORM imaging

% Obtain number of attached probes depending on accuracy of attachment:
numberOfProbes = round(labellingAccuracy*numberOfSurfaceProteins);

% Create empty variables for storing data:
x = [];
y = [];
z = [];
n = []; % number of probes in total


 
% Range to obtain values from =
minvalue = 1.0;
maxvalue = numberOfSurfaceProteins; % number of rows (entries)


for i = 1:numberOfProbes
    
    randomvalue = minvalue+rand(1)*(maxvalue-minvalue);
    
    
    % random value rounded up
    index = round(randomvalue);
    
    % extract entry number 'randomvalue' from the xPoint,
    % yPoint and zPoint array values, and add as probe co-ordinates:
    x = [x; xPoint(index)];
    y = [y; yPoint(index)];
    z = [z; zPoint(index)];

end

n = numberOfProbes;

