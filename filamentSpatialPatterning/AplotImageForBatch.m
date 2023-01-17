function [image, imageArray] = AplotImageForBatch(inputFile, channelColour,xImageSize,yImageSize,scaleIncrease)
% Plot the data from a channel to be analysed.
dataTable = readtable(inputFile);
channelTruthTable = dataTable.Channel == channelColour; %channelColour
figure('position',[10.*scaleIncrease,10.*scaleIncrease,xImageSize.*scaleIncrease,yImageSize.*scaleIncrease])
xlim([0 (50000*scaleIncrease)]);
ylim([0 (80000*scaleIncrease)]);

hold on
%alphaVector = dataTable.Photons/max(dataTable.Photons);
%image = scatter(dataTable.X_nm_(channelTruthTable).*scaleIncrease, dataTable.Y_nm_(channelTruthTable).*scaleIncrease, 0.001, '.k', 'MarkerFaceAlpha', 'flat', 'AlphaData', alphaVector);   %0.01
image = scatter(dataTable.X_nm_(channelTruthTable).*scaleIncrease, dataTable.Y_nm_(channelTruthTable).*scaleIncrease, 0.001, '.k', 'MarkerFaceAlpha', 0.02, 'MarkerEdgeAlpha', 0.02);   %0.01
box off
set(gca,'visible','off')
set(gca,'xtick',[])

%remove the axis for analysis
set(gca,'visible','off')

%Change the m<val> to increase magnification
export_fig('V:/Virus Group/Papers/Vitro Filaments/Figure 4 External proteins/Software/tempFiles/tempImage', '-png', '-m10');
imageArray = print2array(gcf);

end
