clear

inputFolder = 'V:\Andrew\20200930_filamentHAM1\MultiAcquisitionNegative';
plotFigure = false;
numberOfBins = 20;

cd(inputFolder)
imageFiles = dir(strcat(inputFolder, '\\**\\*.tif'));
numberOfImages = length(imageFiles);
lengths = [];

for imageNumber = 1:numberOfImages
    plotImage = imread(strcat(imageFiles(imageNumber).folder(), '\\', imageFiles(imageNumber).name()));
   
    % Need to scale to particular 
    %Skeleton only works off a binary image with the background in black
    adjustedImage = imadjust(plotImage, [0.05 0.20]);     %[0.05 0.20]
    BWPlot = imbinarize(adjustedImage, 0.001);        %0.001 

    %RemoveBackground
    BW2 = bwpropfilt(BWPlot,'Area',[20 999]);
    BW3 = bwmorph(BW2, 'clean');
    BW4 = bwmorph(BW3, 'close');

    skeleton = bwskel(BW4);
    labelledFilaments = bwlabel(skeleton);

    %Don't count Skeletons that have branches
    branchpoints = bwmorph(skeleton, 'branchpoints');
    filamentsWithBranches = unique(labelledFilaments(branchpoints));
    if filamentsWithBranches>0
        labelledFilaments(ismember(labelledFilaments,filamentsWithBranches)) = 0;
    end
    
    %Find the lenths of all filaments in the image
    areas = regionprops(labelledFilaments, 'area');
    lengths = [lengths, areas.Area];
    
    if plotFigure == true
        if max(max(labelledFilaments))>0
            figure
            imshow(adjustedImage)
            xlim([429 856])
            ylim([51 478])
            hold on
            rectangle('Position',[780 450 42.735 7],'FaceColor',[1 1 1], 'EdgeColor', [1 1 1])
            figure
            imshow(BWPlot)
            xlim([429 856])
            ylim([51 478])
            hold on
            rectangle('Position',[780 450 42.735 7],'FaceColor',[1 1 1], 'EdgeColor', [1 1 1])
            figure
            imshow(BW3)
            xlim([429 856])
            ylim([51 478])
            hold on
            rectangle('Position',[780 450 42.735 7],'FaceColor',[1 1 1], 'EdgeColor', [1 1 1])
            %overlay filaments on top of greyscale image
            figure
            imshow(labeloverlay(adjustedImage, labelledFilaments,'Transparency',0));
            xlim([429 856])
            ylim([51 478])
            hold on
            rectangle('Position',[780 450 42.735 7],'FaceColor',[1 1 1], 'EdgeColor', [1 1 1])
            
            %set(gca,'Xlim',[700 1150],'YLim',[200 500])
            %hold on
            %scaleBarSize = 5000/117;
            %rectangle('Position',[1075-(scaleBarSize/2) 450 scaleBarSize 8],'FaceColor',[1 1 1])
            %set(gcf, 'Position',  [10, 10, 2*450, 2*300])
            %set(gca,'visible','off')
        end
    end
    
end

edges = 0:117:3978;
lengths(lengths == 0) = [];
lengths = lengths.*117;
[Y,X] = hist(lengths',edges);
numberOfFilamentsPerFOV = length(lengths')/(numberOfImages*2)
%Y=Y./2;
figure
bar(X,Y,1)
xlim([234 4000])
set(gca,'TickDir','out');
box off
ax=gca;
axes('position',ax.Position,'box','on','ytick',[],'xtick',[],'color','none')
xlim([234 3978])
xlabel('Filament Length (nm)')
ylabel('Frequency')
ylim([0 9700])

%Edit graph vectors to half numbers to allow for 2 colours
