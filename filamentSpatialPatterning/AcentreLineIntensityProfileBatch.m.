%Requires a new version of Matlab for AlphaData command
clear
clc

%Colour of the channels set
green = 0;
red = 1;

%Parameters
colourToBePlotted = red;
scaleIncrease = 0.01; %0.01 gets the whole image on the screen
nmX = 117.0*420.0;
nmY = 117.0*684.0;

%Loop through FOVs
importFolder = 'V:\\Virus Group\\Papers\\Vitro Filaments\\Figure 4 External proteins\\NADataUsed\\red\\';
files = subdir(importFolder);

%Initialise variables for end
averageWidthOfFilaments = []; %create array to hold average widths
allFP1 = zeros(1,1001); %create array to hold fourier transform of filament proteins.  Max 15000nm filament - half length and 2nm sampling

for fileNumber=1:size(files)
    
    %Get file names
    inputFile = files(fileNumber).name;
    
    [inputImage, imageArray] = AplotImageForBatch(inputFile,colourToBePlotted,nmX,nmY,scaleIncrease);
    imageFromFile = imread('V:/Virus Group/Papers/Vitro Filaments/Figure 4 External proteins/Software/tempFiles/tempImage.png');
    greyscaleImage = im2gray(imageFromFile);
    bwImage=imbinarize(greyscaleImage,0.1);   %0.1
    bwImage2=imcomplement(bwImage);

    [imageHeight, imageWidth] = size(greyscaleImage);
    pixelToNmConversion =  nmY/imageHeight;

    %Find centreline/skeleton and perform morphological operations
    se = strel('disk',10);   %10
    closed = imclose(bwImage2,se);
    skeleton = bwskel(closed,'MinBranchLength',50);  %50
    skeleton(bwmorph(skeleton,'branchpoints')) = 0;
    se2 = strel('disk',2);    
    boldSkeleton = imdilate(skeleton, se2);          %For figures, we want bolder skeletons

    %Overlay the centerline with the bwimage
    %{
    figure
    imshow(labeloverlay(greyscaleImage, single(boldSkeleton),'Transparency',0));
    set(gcf,'position',[0,0,nmX,nmY])
    xlim([1600 2850]);
    ylim([6700 7450]);
    %}

    %Label segments
    [labeledImage, numberOfObject] = bwlabel(skeleton);
    [boldLabeledImage, numberOfObject] = bwlabel(boldSkeleton);

    %{
    figure
    imshow(labeloverlay(greyscaleImage, single(boldLabeledImage),'Transparency',0));
    set(gcf,'position',[0,0,nmX,nmY])
    xlim([1600 2850]);
    ylim([6700 7450]);
    %}

    branchMeasurements = regionprops(labeledImage, 'all');
    numberOfBranches = size(branchMeasurements, 1);
    
    if numberOfBranches == 0
        continue
    end

    %Find distance transform
    inverseEdist=bwdist(~bwImage);

    %Find x and y coordinates
    x=[]; y=[]; 
    clear thisBranchPixels; 
    clear branchRDist;
    for k=1:numberOfBranches  
        [row,column]=find(labeledImage==k);
        x{k}=column'; y{k}=row';

        %Find the pixel location
        thisBranchPixels{k} = branchMeasurements(k).PixelIdxList;
        branchRDist{k} = inverseEdist(thisBranchPixels{k});
    end

    %Fit quadratic to points before and after the current point along the 
    % centerline, and look for a perpendicular line 
    pixelsInEachDirection = 10;
    warning('off')
   
    allWidths = zeros(1, length(x));
    
    for j=1:length(x)   % Loop through all the filaments
        lastIndex = length(x{j});
        totalIntensity=zeros((lastIndex),1); %preallocate memory
        widthOfFilament = zeros(lastIndex,1);       %preallocate memory

        if lastIndex>50 && lastIndex<2400       %No particles and nothing too large
            for i=1:lastIndex      % Loop through the pixels in the filament

                index1 = max(1,i - pixelsInEachDirection);
                index2 = min(lastIndex, i + pixelsInEachDirection);

                x1=x{j}(index1:index2);
                y1=y{j}(index1:index2);

                %Fnd the polynomial in this window along the filament
                quadraticAlongFilament = polyfit(x1,y1,2);

                %Slope of line perpendicular to curve
                perpendicularSlope=-1.0./(2.0.*quadraticAlongFilament(1).*x{j}(i)+quadraticAlongFilament(2));

                %The values a radial distance from the centreline are denoted right and
                %left
                rightEndpointXValue=x{j}(i)+branchRDist{j}(i)./sqrt(1.0+perpendicularSlope.^2.0);
                rightEndpointYValue=y{j}(i)+perpendicularSlope.*(rightEndpointXValue-x{j}(i));

                leftEndpointXValue=x{j}(i)-branchRDist{j}(i)./sqrt(1.0+perpendicularSlope.^2.0);
                leftEndpointYValue=y{j}(i)+perpendicularSlope.*(leftEndpointXValue-x{j}(i));

                %warning('on')
                hold off

                %Intensity along filament over line profiles over radii
                Ix=double([leftEndpointXValue rightEndpointXValue]);
                Iy=double([leftEndpointYValue rightEndpointYValue]);
                
                if isnan(Ix(1)) || isnan(Ix(2)) || isnan(Iy(1)) || isnan(Iy(2))
                    continue
                end
                
                q=improfile(greyscaleImage,Ix,Iy);
                totalIntensity(i)=sum(q);  

                %Find the width of the filament at the point being tested
                widthOfFilament(i) = sqrt((leftEndpointXValue - rightEndpointXValue).^2 + (leftEndpointYValue - rightEndpointYValue).^2);

            end

            % Rescale large data points so the don't obscure the data
            totalIntensity(totalIntensity>1000) = 1000;

            xInPixels = [1:lastIndex].*pixelToNmConversion;

            % Plotting
            %{
            figure
            plot(xInPixels,totalIntensity);
            xlabel('Pixels')
            ylabel('Intensity')
            %}
            
            resampledFP1 = A_FFT_Function(xInPixels, totalIntensity);
            
            if isnan(resampledFP1(1))
                continue
            end

            allFP1 = allFP1 + resampledFP1;
            allWidths(j) = mean(widthOfFilament);
        end   
    end

    averageWidthOfFilaments{fileNumber} = allWidths;
    
end

%{
%Final plotting of filament width
averageWidthOfFilaments =  averageWidthOfFilaments.*pixelToNmConversion;
figure
histogram(averageWidthOfFilaments,50);

%Final plotting of Fourier transform
figure
plot([0:0.5:500],allFP1)
xlabel('Frequency/nm^-^1')
ylabel('Normalised intensity')

ylabel('Normalised intensity')
%}
