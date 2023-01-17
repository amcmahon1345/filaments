function xyz2DDStormplot = make2DMeanSStormPlotWithTitle(virion_type, SP1LabelColour, SP2LabelColour, plotProbesAndProteins, numberOfSurfaceProteins, x3,y3,z3, x4,y4,z4)

if plotProbesAndProteins == "true" && numberOfSurfaceProteins > 0
    
    figure
    
    if virion_type == "sphere"
        plot(x3,y3, SP1LabelColour, x4,y4, SP2LabelColour)
        % plot(X,Y) creates a 2-d plot
    
    else
        plot(x3,z3,SP1LabelColour, x4,z4,SP2LabelColour)
    end
    
    set(gca,'Color','k') % black background
    
    daspect([1 1 1]) % needed for scalebar
    scalebar() % default = white
end
