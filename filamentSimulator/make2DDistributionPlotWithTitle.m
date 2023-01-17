function xy2Dplot = make2DDistributionPlotWithTitle(virion_type, SP1LabelColour, SP2LabelColour, plotProbesAndProteins,numberOfSurfaceProteins,x1,y1,z1, x2,y2,z2, x3,y3,z3, x4,y4,z4)

if plotProbesAndProteins == "true" && numberOfSurfaceProteins > 0
    
    figure
    
    if virion_type == "sphere"
        plot(x1,y1,'+', x2,y2,'_m', x3,y3, SP1LabelColour, x4,y4, SP2LabelColour)
    
    else
        plot(x1,z1,'+', x2,z2,'_m', x3,z3,SP1LabelColour, x4,z4,SP2LabelColour)
    end
    
    
    daspect([1 1 1]) % needed for scalebar
    scalebar('Colour', [0 0 0]) % black bar
    
elseif plotProbesAndProteins == "true" % ie if zero surface proteins are present:
        
    plot3(xProbe,yProbe,zProbe,'x')

end
