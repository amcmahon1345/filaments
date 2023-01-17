function [x,y,z] = produceRandomPointsOnTheSurfaceOfAFilamentWithExclusionRadius(numberOfPointsToAdd, exclusionRadiusForPoints, diameterOfFilament, lengthOfFilament, xOtherPoints, yOtherPoints, zOtherPoints, exclusionRadiusForOtherPoints, alternatingProteinLocations, proteinNumber)

% Produce a random number and relative surface areas to see whether a point
% is on the caps or the main body of a filament and then randomly produce a
% point and ensure that it is a minimum distance from other points
    

% Preallocation of zeros for speed and index the point being added
i = 1;
previousNumberOfPoints = size(xOtherPoints);
[x,y,z] = deal(zeros(numberOfPointsToAdd, 1));

xOtherPoints = [xOtherPoints; zeros(numberOfPointsToAdd, 1)];
yOtherPoints = [yOtherPoints; zeros(numberOfPointsToAdd, 1)];
zOtherPoints = [zOtherPoints; zeros(numberOfPointsToAdd, 1)];

exclusionRadiusForOtherPoints = [exclusionRadiusForOtherPoints;
    
    exclusionRadiusForPoints.*ones(numberOfPointsToAdd, 1)];

% til now making list of all other positions, adding zero to all = exclusion list

while x(numberOfPointsToAdd) == 0 && y(numberOfPointsToAdd) == 0 && z(numberOfPointsToAdd) == 0
    
    randomNumber = rand(1);
    
    % placing proteins on cylinder/ caps
    % probabiliy of putting somehting on cylinder = proportion to SA of
    % cylinder, equally for caps. (caps= same as sphere) for sphere =
    % 4pir^2. Then surface area of cylinder = 2pir*length. therefore
    % 2r/length. Therefore if rand(1) (between 0 and 1) is smaller than the
    % ratio then thats the probabilty of being on cylinder (body). Ie if
    % diameter is smaller therefore no body and so on cap. basically ratio
    % of surface areas to determine where it is.

    if randomNumber < ((lengthOfFilament-diameterOfFilament)/lengthOfFilament)
        
        %For Normal
        if alternatingProteinLocations == 0
            [x(i),y(i),z(i)] = produceRandomPointsOnTheSurfaceOfTheFilamentBodyWithExclusionRadius(exclusionRadiusForPoints, diameterOfFilament, lengthOfFilament, xOtherPoints, yOtherPoints, zOtherPoints, exclusionRadiusForOtherPoints);
        
        %For alternating filament
        else
            if proteinNumber == 1
                [x(i),y(i),z(i)] = produceSinusoidRandomPointsOnTheSurfaceOfTheFilamentBodyWith(exclusionRadiusForPoints, diameterOfFilament, lengthOfFilament, xOtherPoints, yOtherPoints, zOtherPoints, exclusionRadiusForOtherPoints);
            else
                [x(i),y(i),z(i)] = produceCosineRandomPointsOnTheSurfaceOfTheFilamentBodyWith(exclusionRadiusForPoints, diameterOfFilament, lengthOfFilament, xOtherPoints, yOtherPoints, zOtherPoints, exclusionRadiusForOtherPoints);
            end
        end
    else
        [x(i),y(i),z(i)] = produceRandomPointsOnTheSurfaceOfTheCapsWithExclusionRadius(exclusionRadiusForPoints, diameterOfFilament, lengthOfFilament, xOtherPoints, yOtherPoints, zOtherPoints, exclusionRadiusForOtherPoints);
    end
    
    xOtherPoints(i + previousNumberOfPoints) = x(i);
    yOtherPoints(i + previousNumberOfPoints) = y(i);
    zOtherPoints(i + previousNumberOfPoints) = z(i);
    
    i=i+1;
end
