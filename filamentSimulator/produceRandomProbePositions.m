function [x,y,z] = produceRandomProbePositions(n, numberOfProbes, muX, muY, muZ, FWHM);
% Obtaining random normal distribution of probe locations on
% surface protein for STORM image plot


% INPUT:
% numberOfProbes would = numberOfHA/numberOfNA
% muX/muY = mean = the xHAprobe pos/yHAprobe pos etc
% sigma = standard deviation = FWHM/2.355
% n = number of points wanted (this would depend on number of flashes - can
% be obtained from previous real data as discussed with Andrew)

sigma = FWHM/2.355;

% Create empty variables for storing data:
x = [];
y = [];
z = [];

for i = 1:numberOfProbes % forloop going through each probe
    
    % Using randn() to get the random positions
    total1 = randn(n, 1); % so total points = n
    total2 = randn(n, 1);
    total3 = randn(n, 1);% to get different random numbers for the x and 
                % y to allow sampling of whole 2d space otherwise just a
                % straight line
                
    
    
    for j = 1:n % for that probe want to create new postions 
                % therefore forloop going through
                    
        % extract value from 'total' variable
        index1 = total1(j);
        index2 = total2(j);
        index3 = total3(j);
    
        % each position of the 'total' variable multiplied by both mus
        valueX = index1*sigma + muX(i); 
        valueY = index2*sigma + muY(i);
        valueZ = index3*sigma + muZ(i);
    
        % gives new pos (valueX/Y) which is added to the empty variable
        x = [x; valueX]; 
        y = [y; valueY];
        z = [z; valueZ];
    
    end % end of one probe variations
    
    
end % gone through all probes and made variations
