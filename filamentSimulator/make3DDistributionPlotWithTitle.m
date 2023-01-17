function make3DDistributionPlotWithTitle(plotProbesAndProteins,SP1LabelColour, SP2LabelColour,numberOfSurfaceProteins,x1,y1,z1, x2,y2,z2, x3,y3,z3, x4,y4,z4)

    if plotProbesAndProteins == "true" && numberOfSurfaceProteins > 0

        figure
        plot3(x1,y1,z1,"+", x2,y2,z2,"_m", x3,y3,z3, SP1LabelColour, x4,y4,z4, SP2LabelColour)
        daspect([1 1 1]) % needed for scalebar

    elseif plotProbesAndProteins == "true" % ie if zero surface proteins are present:

        figure
        plot3(x1,y1,z1,'x', x2,y2,z2,'o')
        daspect([1 1 1]) % needed for scalebar

    end
