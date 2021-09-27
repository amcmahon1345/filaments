function [x,y,z] = produceRandomPointsOnTheSurfaceOfTheFilamentBodyWithExclusionRadius(exclusionRadiusForPoint, diameterOfFilament, lengthOfFilament, xOtherPoints, yOtherPoints, zOtherPoints, exclusionRadiusForOtherPoints);
% Producing random points on the surface of the body of a filament

pointPlaced = false; % flag to say if greater than exclusion radius its been added to list and model.

while pointPlaced == false
    randomNumbers = rand(1, 2);% making 2 random numbers  
    tempZ = randomNumbers(:,1)*(lengthOfFilament-diameterOfFilament); % first for distance along the cylinder (z - pos)
    theta = randomNumbers(:,2).*pi.*2; % this is z-pos
    
    % position of a random point:
    [tempX,tempY,tempZ] = pol2cart(theta, diameterOfFilament/2, tempZ); 
    
    % pol2cart = a function.
    % theta, radial and z co-ord = input.
    % using this function we obtain our "x" and "y" temporary positions.
    
    
    % distance to all other points placed on cylinder (body):
    distance = sqrt((xOtherPoints-tempX).^2 + (yOtherPoints-tempY).^2 + (zOtherPoints-tempX).^2); 
    
    if distance > exclusionRadiusForOtherPoints+exclusionRadiusForPoint 
        % if outside exclusion radius then the temp x, temp y, and temp z 
        % becomes actual x, y, and z positions.
        % causing pointplaced to be true and ending the function.
        
        x = tempX;
        y = tempY;
        z = tempZ;  
        pointPlaced = true;
    end
    
end