function [scaleFourierTransform, intensityFourierTransform] = make2DDStormPlotWithTitleAndFourier(virion_type, SP1LabelColour, SP2LabelColour, plotProbesAndProteins, numberOfSurfaceProteins, numberOfProbes1, x3,y3,z3, numberOfProbes2, x4,y4,z4, n, fwhm, scaleFourierTransform, intensityFourierTransform)

if plotProbesAndProteins == "true" && numberOfSurfaceProteins > 0
    
    [x3r, y3r, z3r] = produceRandomProbePositions(n, numberOfProbes1, x3, y3, z3, fwhm);
    [x4r, y4r, z4r] = produceRandomProbePositions(n, numberOfProbes2, x4, y4, z4, fwhm);
    
    %figure
    
    if virion_type == "sphere"
        plot(x3r,y3r, SP1LabelColour, x4r,y4r, SP2LabelColour)
    
    else
        %plot(x3r,z3r, SP1LabelColour, x4r,z4r, SP2LabelColour)
    end
 
    set(gca,'Color','k') % black background
    daspect([1 1 1]) % needed for scalebar
    scalebar() % default is white
    
    %Code to make histograms for fourier transforms
    datapoints = 2000; 
    [intensityOverProbeLength, ~] = hist(z3r,datapoints);
    maxZ = max(z3r);
    minZ = min(z3r);
    zCoordinates = linspace(minZ, maxZ, datapoints);
    [f, P1] = A_FFT_Function(zCoordinates, intensityOverProbeLength, 1/(zCoordinates(2)-zCoordinates(1)));
    %P1 = interp1(f, P1, [min(f):(f(3)-f(2))/10:max(f)], 'cubic');
    %f = [min(f):(f(3)-f(2))/10:max(f)];
    intensityFourierTransform = intensityFourierTransform + P1;
    scaleFourierTransform = f;
    
end