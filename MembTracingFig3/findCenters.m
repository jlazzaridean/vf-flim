%This function accepts a mask of traced membranes as well as
%minimum and maximum areas for the cell centers to find. Good values of
%these parameters seem to be 10000 and ~2 respectively (to remove the
%background and debris) with the zoom 0.7 images. The min size threshold
%needs to be a lot smaller than the minArea for the main
%traceMembranes function, as the cells are a lot more "filled in" when
%they are processed at this stage.

%This function dilates the mask to close up as many small gaps in the
%membrane as possible. It then looks for the number of non-connected cell
%centers (inverse of the mask). The counting is fairly robust, but it
%undercounts where the membrane is discontinuous and the ROI is open to the
%background. 

%The number of
%non-connected groups in the correct size range is returned as the number
%of cells.

%Last edits: Julia Lazzari-Dean, September 3, 2018

function [nCells] = findCenters(mask,maxArea,minArea)

%dilate the mask with a diamond structuring element to 'close up' some of
%the gaps where membranes go below threshold
dilMask = imdilate(mask,strel('diamond',2));

%invert to get cell centers
centers = imcomplement(dilMask);

%map the connectivity and remove both debris and the background
cc = bwconncomp(centers);
%find the properties of the regions that were identified
stats = regionprops(cc);
%remove the background (bigger than the max area)
tooLarge = [stats.Area]>maxArea;
%remove the debris (smaller than the min area)
tooSmall = [stats.Area]<minArea;

%subtract the areas that are too small or too large from the total count.
nCells = cc.NumObjects - sum(tooLarge) - sum(tooSmall);

end