% Make a virus with randomly distributed proteins and add a number of probes

clear   % remove items from workspace
clc     % clear command window

% Input virion parameters:
virion_type = "filament"; % type of virion either "sphere" or "filament"
diameterOfVirus = 80;  % diameter size in "nm"
lengthOfVirus = 500;    % If the virus is spherical, this is also the virus diameter
numberOfTypesOfSurfaceProteins = 1; % max = 2 surface proteins
alternatingProteinLocations = 0;    % 1 is alternating with sinusoid- probe 2 out of phase, 0 is non-alternating
numberOfSurfaceProteins = 500;   %375
exclusionRadiusOfSurfaceProteins = 0;   % 3.5 is the size of the radius of the head of nCov-19

% If two surface proteins,split surface proteins(SP):
minSP1toSP2ratio = 6.02; % ie NA:HA = 1:6.02 
maxSP1toSP2ratio = 7.63; % ie NA:HA 7.63

% Input probe/labelling parameters:
labellingEfficiency = 1; % ie. 0.1 = 10% accuracy
SP1LabelColour = ".g"; % '.' = a dot in the plot and letter for colour
SP2LabelColour = ".r"; % g for green and r for red
exclusionRadiusOfProbes = 0;

% Input imaging parameters:
numberOfSimulations = 1;
plotProbesAndProteins = "true";
flashes = 50; % relying on duration of imaging to determine input value
FWHM = 2.355 * 7.4; % aka optical resolution impacting sigma of flashes (nm)

% Initialised variables
intensityFourierTransform = zeros(1,1001);  
scaleFourierTransform = zeros(1,1001);   

% Compute:
for i = 1:numberOfSimulations       % Loop over as many iterations of simulations required
    
    % Only need to place surface proteins if there are any:
    if numberOfSurfaceProteins ~= 0   % '~=' refers to 'not equal to'
        
        % Obtain exact number of surface proteins
        if numberOfTypesOfSurfaceProteins == 1
            
            numberOfSurfaceProtein1 = numberOfSurfaceProteins;
            
            % Surface protein co-ordinates:
            [xSP1, ySP1, zSP1] = produceRandomPointsOnTheSurfaceOfAFilamentWithExclusionRadius(numberOfSurfaceProtein1, exclusionRadiusOfSurfaceProteins, diameterOfVirus, lengthOfVirus, [], [], [], [], alternatingProteinLocations, 1); % empty brackets for allowing input of other probe locations, so that

            % Probe/Label co-ordinates: 
            [xSP1Probe, ySP1Probe, zSP1Probe, SP1Probes] = obtainProbePositions(labellingEfficiency, numberOfSurfaceProtein1, xSP1, ySP1, zSP1);   
        
            [xSP2, ySP2, zSP2, xSP2Probe, ySP2Probe, zSP2Probe, SP2Probes] = deal(0);
            
            
            % Plotting to check distribution of proteins and/or probes:
            %make3DDistributionPlotWithTitle(plotProbesAndProteins, SP1LabelColour, SP2LabelColour, numberOfSurfaceProteins, xSP1, ySP1, zSP1, xSP2, ySP2, zSP2, xSP1Probe, ySP1Probe, zSP1Probe, xSP2Probe, ySP2Probe, zSP2Probe)
            % 'plus' in plot = SP1
            % red point = centre
            
            % 2-D Distribution Plot:
            %make2DDistributionPlotWithTitle(virion_type, SP1LabelColour, SP2LabelColour, plotProbesAndProteins, numberOfSurfaceProteins, xSP1, ySP1, zSP1, xSP2, ySP2, zSP2, xSP1Probe, ySP1Probe, zSP1Probe, xSP2Probe, ySP2Probe, zSP2Probe)
        
            
            % 2-D DStorm Simulation Plots:
            % a) showing exact location with one flash:
            %make2DMeanDStormPlotWithTitle(virion_type, SP1LabelColour, SP2LabelColour, plotProbesAndProteins, numberOfSurfaceProteins, xSP1Probe, ySP1Probe, zSP1Probe, xSP2Probe, ySP2Probe, zSP2Probe);
        
            % b) showing random normal distribution flashes:
            [scaleFourierTransform, intensityFourierTransform] = make2DDStormPlotWithTitleAndFourier(virion_type, SP1LabelColour, SP2LabelColour, plotProbesAndProteins, numberOfSurfaceProteins, SP1Probes, xSP1Probe, ySP1Probe, zSP1Probe, SP2Probes, xSP2Probe, ySP2Probe, zSP2Probe, flashes, FWHM, scaleFourierTransform, intensityFourierTransform);
            
            
        else % two surface proteins present = 
            [numberOfSurfaceProtein1, numberOfSurfaceProtein2] = produceNumberOfSurfaceProtein(minSP1toSP2ratio, maxSP1toSP2ratio, numberOfSurfaceProteins);
        
            % Surface protein co-ordinates:
            [xSP1, ySP1, zSP1] = produceRandomPointsOnTheSurfaceOfAFilamentWithExclusionRadius(numberOfSurfaceProtein1, exclusionRadiusOfSurfaceProteins, diameterOfVirus, lengthOfVirus, [], [], [], [], alternatingProteinLocations, 1); % empty brackets for allowing input of other probe locations, so that 
            [xSP2, ySP2, zSP2] = produceRandomPointsOnTheSurfaceOfAFilamentWithExclusionRadius(numberOfSurfaceProtein2, exclusionRadiusOfSurfaceProteins, diameterOfVirus, lengthOfVirus, [], [], [], [], alternatingProteinLocations, 2);
           
            % Probe/Label co-ordinates: 
            [xSP1Probe, ySP1Probe, zSP1Probe, SP1Probes] = obtainProbePositions(labellingEfficiency, numberOfSurfaceProtein1, xSP1, ySP1, zSP1);
            [xSP2Probe, ySP2Probe, zSP2Probe, SP2Probes] = obtainProbePositions(labellingEfficiency, numberOfSurfaceProtein2, xSP2, ySP2, zSP2);
        
        
            % Plotting to check distribution of proteins and/or probes:
            make3DDistributionPlotWithTitle(plotProbesAndProteins, SP1LabelColour, SP2LabelColour, numberOfSurfaceProteins, xSP1, ySP1, zSP1, xSP2, ySP2, zSP2, xSP1Probe, ySP1Probe, zSP1Probe, xSP2Probe, ySP2Probe, zSP2Probe)
            % 'plus' in plot = SP1 and 'minus' = SP2
            
            % 2-D Distribution Plot:
            make2DDistributionPlotWithTitle(virion_type, SP1LabelColour, SP2LabelColour, plotProbesAndProteins, numberOfSurfaceProteins, xSP1, ySP1, zSP1, xSP2, ySP2, zSP2, xSP1Probe, ySP1Probe, zSP1Probe, xSP2Probe, ySP2Probe, zSP2Probe)
        
            
            % 2-D DStorm Simulation Plots:
            % a) showing exact location with one flash:
            make2DMeanDStormPlotWithTitle(virion_type, SP1LabelColour, SP2LabelColour, plotProbesAndProteins, numberOfSurfaceProteins, xSP1Probe, ySP1Probe,zSP1Probe, xSP2Probe,ySP2Probe,zSP2Probe);
        
            % b) showing random normal distribution flashes:
            [scaleFourierTransform, intensityFourierTransform] = make2DDStormPlotWithTitleAndFourier(virion_type, SP1LabelColour, SP2LabelColour, plotProbesAndProteins, numberOfSurfaceProteins, SP1Probes, xSP1Probe, ySP1Probe,zSP1Probe, SP2Probes, xSP2Probe,ySP2Probe, zSP2Probe, flashes, FWHM, scaleFourierTransform, intensityFourierTransform);
        end
    end 
end 

figure
plot(scaleFourierTransform, intensityFourierTransform)