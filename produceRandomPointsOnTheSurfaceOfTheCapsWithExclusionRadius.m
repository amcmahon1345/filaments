function [x,y,z] = produceRandomPointsOnTheSurfaceOfTheCapsWithExclusionRadius(exclusionRadiusForPoint, diameterOfFilament, lengthOfFilament, xOtherPoints, yOtherPoints, zOtherPoints, exclusionRadiusForOtherPoints)
% Producing random points on the surface of a sphere with an exclusion
% distance around each probe


pointPlaced = false;

while pointPlaced == false
    randomNumbers = rand(1, 2);
    
    % now looking on surface of sphere
    % find pos on one sphere, split in half and by the length of filament
    % body you move it across
    
    % spherical polar co-ords:
    
    phi = acos(2.*randomNumbers(:,1) -1) - (pi/2);  % acos to get poles to have same density. 
                                                    % pi/2 as acos returns positive values
                                                    % obtaining pos with
                                                    % equal probabilty
                                                    % across hemisphere
                                                    % rather than at poles.
                                                    
                                                    % search acos fn and
                                                    % also goodole placing
                                                    % random points on
                                                    % sphere.
                                                    
                                                   
    theta = randomNumbers(:,2).*pi.*2; % same as on body
    
    radius = (diameterOfFilament/2);
    
    [tempX,tempY,tempZ] = sph2cart(theta, phi, radius); % position of a random point
    
    % sph2cart = fn used to obtain 
    
    if tempZ > 0 % this if statement is the splitting of the sphere accroding to the middle
        tempZ = tempZ+lengthOfFilament-diameterOfFilament;
    end
    
    % distance to other points placed: basically exclusion radius)
    distance = sqrt((xOtherPoints-tempX).^2 + (yOtherPoints-tempY).^2 + (zOtherPoints-tempX).^2); 
    
    if distance > exclusionRadiusForOtherPoints+exclusionRadiusForPoint
        x = tempX;
        y = tempY;
        z = tempZ;
        pointPlaced = true;
    end
    
end