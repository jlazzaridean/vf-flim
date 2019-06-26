%This is the main function for processing of exported (.asc) lifetime data
%from SPCImage. Written by Julia Lazzari-Dean, August 2018. Dependencies
%are the following matlab functions: calculateStats.m, traceMembranes3.m,
%struct2PandasCSV.m, findCenters.m matchAndCategorize.m, maskImage_otsu.m,
%findROIsToMerge.m, ccToMask.m

%Overall, this function finds cells in photon count
%images using Otsu's method of thresholding. It then identifies isolated
%groups of cells using Matlab's built-in connectivity assessment
%(bwconncomp). From these connectivity maps, it generates regions of
%interest (ROIs) for analysis of resting membrane potential in the
%corresponding lifetime image. The user is asked to approve ROIs generated for
%every image (rainbow color scheme delineates cell groups) by
%pressing "Enter" when the image appears. If any debris hasn't been
%filtered out (i.e. if debris is color-coded), it can be removed by a left
%click on the offending ROI before the image is approved with "Enter."

%The script then applies these photon count image ROIs to the lifetime
%image. The mean value of lifetime in each ROI is tabulated for each image,
%coverslip, and category (entered in "metadata", see below). These values
%are saved to .mat files in an output directory, and they are also exported
%to .csv for input into Pandas dataframes (python) or other convenient
%plotting software.

%The script generates two box plots, one of the data broken down by
%coverslip, and the other of the data broken down by category. These plots
%are also saved to the output directory.

%The input parameters to the function are..
%1.the size of the smallest object to include (in pixel units). 
%With zoom 1 data on the 40x objective, 50 or 70 work well empirically.
%2. metadata, which should be a cell array. The first column should be a
%numerical, unique coverslip ID. The second column should be a string
%describing the category (i.e. 'highK' or 'BIB'). Categories may be
%repeated.
%3. The name of an output folder that will be generated in the directory
%with the images. The output files are put in a separate folder so that the
%same data may be analyzed multiple times without "trampling" the values
%(i.e. just put in a new folder name).
%4. A base name for output files (graphs, .csv, etc.). Enter as a character
%vector in single quotations and do not include a file extension.
%5. A string corresponding to the type of cell that this is (helps with the
%thresholding algorithms - some of the levels used are hardcoded in by cell
%type to help it do better with A431s vs. HEK293T).

%Slightly finicky aspects of this code:
%1. The matching of lifetime and photon count .asc files is done by what is
%sequentially in the selected directory. The script will check that the
%filenames match, but if sequential filenames for some reason are not real
%pairs, it will return errors and not proceed.
%2. All files in the selected directory should have corresponding entries
%in the metadata input. The code is not fully tested for tolerance to
%messing up the input format.
%3. All image files must be in the format "YYYY-MM-DD_XX-XX_ for proper
%processing (where XX-XX is the coverslip number). X-X, X-XX, XX-X (i.e.
%changing the presence or absence of leading zeros) should be ok.

%Written by Julia Lazzari-Dean.
%Last edits 9/3/18

function [] = batchTraceMembranes2(smallestObject,metadata,outputFolderName,outputFileName,cellType)

%prompt the user for the relevant directory and go to that directory
matlab_dir = pwd;
directory = uigetdir(matlab_dir,'Please select folder containing .asc photon and tm files.');
try
    cd(directory);
catch
    disp('No directory selected. Exiting.')
    return
end
%get the correct file separator
fsep = filesep;

%get all .asc files in this directory
filesList = dir('*.asc');
%make a directory to hold the output files and prevent over-writing
mkdir(outputFolderName);
%go back to the matlab directory so dependencies function normally
cd(matlab_dir);

%make a structure to hold the appropriately matched files and metadata, as
%well as the later analysis
coverslipResults = struct('coverslipID',0,'fileList',"",'meanLifetimes',0,'category',"",'nCells',0,'imageID',0);
categorizedResults = struct('category',"",'fileList',"",'meanLifetimes',0);

%instantiate the results structures based on the metadata
%find the categories listed in the metadata cell array
categories = strings(size(metadata,1),1);
for i = 1:size(metadata,1)
    categories(i,1) = metadata{i,2};
end
%only take the unique categories for instantiating the struct
uniqueCategories = unique(categories);
for i = 1:size(uniqueCategories,1)
    categorizedResults(i).category = uniqueCategories(i,1);
    categorizedResults(i).fileList = "";
end

%instantiate the coverslip results structure based on the coverslip IDs
for i = 1:size(metadata,1)
    coverslipResults(i).coverslipID = metadata{i,1};
    coverslipResults(i).fileList = "";
end

%sort the input files into structs based on filename, assuming there is one
%photons and one color coded value image for each
photonList = strings(length(filesList)/2,1);
lifetimeList = strings(length(filesList)/2,1);
photonCount = 1;
lifetimeCount = 1;
%separate files into those containing "photons" and those containing "color
%coded value". If a file does not contain one of these handles, throw an
%error.
for i = 1:length(filesList)
    if (~isempty(strfind(filesList(i,1).name,'_photons.asc')))
        photonList(photonCount,1) = filesList(i,1).name;
        photonCount = photonCount + 1;
    elseif (~isempty(strfind(filesList(i,1).name,'_color coded value.asc')))
        lifetimeList(lifetimeCount,1) = filesList(i,1).name;
        lifetimeCount = lifetimeCount + 1;
    else
        error('Improperly formatted file.');        
    end
end

%make a structure for the lifetime-photon image pairs
imgList = struct('photonName',0,'lifetimeName',0,'baseName',0,'directory',0,'coverslipID',0,'imageID',0,'category',"");
%check that the lifetime and photon count images are 'matched' in number
if (photonCount ~= lifetimeCount)
    error('Unequal numbers of lifetime and photon count images.');
end
numImages = photonCount-1;

%match up the photon count and the flim image and categorize them
for i = 1:numImages
    %find the image names for the photon and lifetime images at position i
    imgList(i).photonName = photonList(i);
    imgList(i).lifetimeName = lifetimeList(i);
    imgList(i).directory = directory;
    
    %call the function that parses the filenames and puts them in
    %categories
    [imgList(i).baseName, imgList(i).category, imgList(i).coverslipID, imgList(i).imageID] = matchAndCategorize(imgList(i).photonName,imgList(i).lifetimeName,metadata);
end

%print out a table with a list of all the grouped files
%removing the parts of imageList that will look redundant for easy viewing
T = struct2table(imgList);
T = removevars(T,{'photonName','lifetimeName','directory'});
disp(T)

%ask the user to approve the groupings
answer = "";
while(~strcmp(answer,'Y') && ~strcmp(answer,'N') && ~strcmp(answer,'n') && ~strcmp(answer,'y'))
    prompt = 'Are these categories what you wanted? [Y/N] ';
    answer = input(prompt,'s');
end  
%process the user input and either exit or continue.
if(strcmp(answer, 'Y') || strcmp(answer,'y'))
    disp('Beginning image analysis with provided categories.');
elseif (strcmp(answer, 'N') || strcmp(answer,'n'))
    error('Categories rejected.');
end        

%open the files, trace them, and save/categorize the results
for i = 1:numImages
    %load the paired image files
    photons = load([char(imgList(i).directory),fsep,char(imgList(i).photonName)]);
    flim = load([char(imgList(i).directory),fsep,char(imgList(i).lifetimeName)]);
    
    %call the tracing script - keep running it until the user approves the
    %ROIs
    goodROIs = 0;
    while(~goodROIs)
        [results, goodROIs] = traceMembranes3(photons,flim,cellType,smallestObject);
    end
    
    %save the results and the file names together in a .mat file in the
    %parent directory for that image
    outputName = [char(imgList(i).directory),fsep,outputFolderName,fsep,char(imgList(i).baseName),'_traced'];
    filePair = imgList(i);
    save(outputName, 'results');
    save(outputName, 'filePair', '-append');

    %add the results to the appropriate coverslip's list in the structure
    %find the correct coverslip in the overall structure
    for j = 1:size(coverslipResults,2)
        if(coverslipResults(j).coverslipID == imgList(i).coverslipID)
            if(coverslipResults(j).fileList == "")
                %start the lists if they doesn't exist yet
                coverslipResults(j).fileList(1,1) = imgList(i).baseName;
                %also add the mean lifetime values - need to iterate
                %through the struct and record a value from each ROI
                for k = 1:size(results,2)
                    if(coverslipResults(j).meanLifetimes == 0)
                        coverslipResults(j).meanLifetimes(1,1) = results(1,k).mean; 
                        coverslipResults(j).nCells(1,1) = results(1,k).nCells;
                        coverslipResults(j).imageID(1,1) = imgList(i).imageID;
                    else
                        coverslipResults(j).meanLifetimes(end+1,1) = results(1,k).mean;
                        coverslipResults(j).nCells(end+1,1) = results(1,k).nCells;
                        coverslipResults(j).imageID(end+1,1) = imgList(i).imageID;
                    end
                end
            else
                %if the lists do exist, add to the end
                coverslipResults(j).fileList(end+1,1) = imgList(i).baseName;
                %also add the mean lifetime values - need to iterate
                %through the struct and record a value from each ROI
                for k = 1:size(results,2)
                    coverslipResults(j).meanLifetimes(end+1,1) = results(1,k).mean;
                    coverslipResults(j).nCells(end+1,1) = results(1,k).nCells;
                    coverslipResults(j).imageID(end+1,1) = imgList(i).imageID;
                end
            end
        end
    end
    
    %close any open figures          
    close all;
end

%put each coverslip's results into the appropriate category based on the
%metadata, and add the category to the coverslip struct
for i = 1:size(coverslipResults,2)
    %find the correct category for the coverslip from the metadata
    category = "";
    for j = 1:size(metadata,1)
        %iterate across the metadata for each coverslip until the
        %coverslipID matches
        if(coverslipResults(i).coverslipID == metadata{j,1})
            %save the string corresponding to this category
            category = metadata{j,2};
        end
    end
    
    %add the category to the coverslip struct
    coverslipResults(i).category = category;
    
    %now look for this string in the categorizedResults struct and add the
    %results in as appropriate
    for j = 1:size(categorizedResults,2)
        if(strcmp(category,categorizedResults(j).category))
            if(categorizedResults(j).fileList == "")
                %instantiate the file list and data list if it hasn't already
                %been done
                categorizedResults(j).fileList = coverslipResults(i).fileList;
                categorizedResults(j).meanLifetimes = coverslipResults(i).meanLifetimes;
            else
                %if the file list and data list already exist, add to them
                %first need to determine how many items there are to
                %correctly index into the existing datastructure.
                nFiles = numel(coverslipResults(i).fileList);
                nMeasurements = numel(coverslipResults(i).meanLifetimes);
                categorizedResults(j).fileList((end+1):(end+nFiles),1) = coverslipResults(i).fileList;
                categorizedResults(j).meanLifetimes((end+1):(end+nMeasurements),1) = coverslipResults(i).meanLifetimes;
            end
        end
    end
end

%pull the data to plot out of the struct with the results
plotCoverslip = zeros(1,size(coverslipResults,2));
Xt = zeros(1,size(coverslipResults,2));
for i = 1:size(coverslipResults,2)
    numResults = numel(coverslipResults(i).meanLifetimes);
    plotCoverslip(1:numResults,i) = coverslipResults(i).meanLifetimes;
    Xt(1,i) = coverslipResults(i).coverslipID;
end
%remove zero values
plotCoverslip(plotCoverslip == 0) = NaN;

%pull the data to plot out of the struct with the results
plotCategory = zeros(1,size(categorizedResults,2));
Xt_cat = strings(1,size(categorizedResults,2));
for i = 1:size(categorizedResults,2)
    numResults = numel(categorizedResults(i).meanLifetimes);
    plotCategory(1:numResults,i) = categorizedResults(i).meanLifetimes;
    Xt_cat(1,i) = categorizedResults(i).category;
end
%remove zero values
plotCategory(plotCategory == 0) = NaN;

%calculate some statistics on each data grouping - do this first so you can
%print the sample sizes on the graph
statsByCoverslip = calculateStats(plotCoverslip);
statsByCategory = calculateStats(plotCategory);
%add a better label to the stats by category data table
for i = 1:size(statsByCategory,2)    
    statsByCategory(i).category = Xt_cat(i);
end
%display these values in a table
disp(struct2table(statsByCoverslip));
disp(struct2table(statsByCategory));

%graph the results by coverslip
f1 = figure;
figList = f1;
ax1 = gca;
axList = ax1;

%make the boxplot by coverslip
bp = boxplot(plotCoverslip);
xlabel('Coverslip ID');
ylabel('Lifetime (ps)');
title([outputFileName ' - by Coverslip'],'interpreter','none'); %this name-value pair disables underscores as subscripts
xticklabels(Xt);
ylim([1000 2000]);

%print the sample sizes for each coverslip on the graph
for i = 1:size(coverslipResults,2)
    szTxt = ['n = ' num2str(statsByCoverslip(i).count)];
    text(i,1100,szTxt,'HorizontalAlignment','center');
end

%graph the results by category
f2 = figure;
figList(end+1) = f2;
ax2 = gca;
axList(end+1) = ax2;

%make the boxplot by category
bp2 = boxplot(plotCategory);
xlabel('Category');
ylabel('Lifetime (ps)');
title([outputFileName ' - by Category'],'interpreter','none'); %this name-value pair disables underscores as subscripts
xticklabels(Xt_cat);
ylim([1000 2000]);

%print the sample sizes for each condition on the graph
for i = 1:size(categorizedResults,2)
    szTxt = ['n = ' num2str(statsByCategory(i).count)];
    text(i,1100,szTxt,'HorizontalAlignment','center');
end


%set the aspect ratio based on the number of categories/coverslips
if(size(coverslipResults,2) > 4)
    pbaspect(ax1,[1.5 1 1]);
else
    pbaspect(ax1,[1 1 1]);
end
if(size(categorizedResults,2) > 4)
    pbaspect(ax2,[1.5 1 1]);
else
    pbaspect(ax2,[1 1 1]);
end

%clean up the figures
for i=1:size(figList,2)
    set(figList(i),'PaperUnits','inches');
    set(figList(i),'PaperSize', [7 6]);
    set(figList(i),'PaperPosition',[1 0.5 6.5 5.5]);
    set(figList(i),'PaperPositionMode','Manual');
    set(figList(i),'color','w');
end

for i=1:size(axList,2)
    set(get(axList(i),'xlabel'),'FontSize', 14,'Color',[0 0 0]);
    set(get(axList(i),'ylabel'),'FontSize', 14,'Color',[0 0 0]);
    set(get(axList(i),'title'),'FontSize', 16,'Color',[0 0 0]);
    set(axList(i),'XColor',[0 0 0],'YColor',[0 0 0],'ZColor',[0 0 0]);
    set(axList(i),'TickDir','in')
    set(axList(i),'box','off');
    set(axList(i),'LineWidth',0.75);
    set(axList(i),'FontSize',14);
end

%save the results as .mat files
overallOutputName = strcat(directory,fsep,outputFolderName,fsep,outputFileName);
save(overallOutputName, 'categorizedResults');
save(overallOutputName, 'coverslipResults', '-append');
save(overallOutputName, 'metadata', '-append');
save(overallOutputName, 'smallestObject', '-append');
save(overallOutputName, 'statsByCoverslip', '-append');
save(overallOutputName, 'statsByCategory', '-append');

%save the plots!
f1Name = strcat(overallOutputName,'_byCoverslip');
saveas(f1,[char(f1Name),'.pdf']);
f2Name = strcat(overallOutputName,'_byCat');
saveas(f2,[char(f2Name),'.pdf']);

%output the data from the coverslip and category datasets to .csv for
%easier opening in python
outPath = strcat(directory,fsep,outputFolderName);
struct2PandasCSV(coverslipResults,["meanLifetimes","nCells","imageID"],["coverslipID","category"],char(outPath),char(outputFileName));

end