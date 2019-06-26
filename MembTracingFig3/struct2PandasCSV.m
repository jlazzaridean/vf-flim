%This function converts a Matlab structure (dataStructure) into a expanded
%.csv file ready to import into Python's Pandas. The general structure of
%this file will be each row is a measurement with its corresponding
%modifiers (e.g. conditions) that can be read as categorical data and
%easily manipulated.

%The input structure must contain at least one field ("longFields") that is an array of values.
%The must be other "modifier fields" that need to be associated with these
%values. For example, longFields might be lifetime values, and the
%modifier fields might be the coverslip or condition that the value was
%taken from. For each entry in the struct, there can be multiple values in
%the "longest field," but there must be a single unique value in the
%modifier fields. 

%If multiple types are passed to longest field, they must be the same
%lenght (i.e. take the meanLifetime and the number of cells for each
%ROI). They must be passed in as an array of strings.

%The script iterates over the structure and pulls out the values from the
%"longest field" array into a column of an output variable. The other
%columns in the output table are the "modifier fields", or metadata from
%the struct that should be exported with the data. The script can handle
%multiple modifier fields (to be passed in a 1XN an array of strings).

%The resulting reshaped data structure is written to the folder outputPath
%under the file name outputName (pass both in as character arrays; the
%.csv extension is added by the program - don't add it).

%It was originally intended to be used with 'batchTraceMembranes' so that
%the outputs of that program can be written to .csv files for easy input
%into Python for plotting and analysis to avoid copy-paste time and errors.

%Sample input line:
%struct2PandasCSV(testerCat,"meanLifetimes",["category","otherModifier"],'C:\Users\Julia\PhD\ImagingData\Deckard\2018-08-01_highKregularBIB_293t\fluorescein_zscore','testMe')
%where meanLifetimes, category, and otherModifier are all fields in a
%testerCat structure.

%Written by Julia Lazzari-Dean
%Last edited August 21,2018

function [] = struct2PandasCSV(dataStructure,longFields,modifierFields,outputPath,outputName)

%determine how many parameters are being processed
nMod = size(modifierFields,2);
nLong = size(longFields,2);

%set up a destination cell array with the correct number of columns
nCols = nMod + nLong;
outputData = cell(1,nCols);

%put headers on the csv file for the longest field modifiers
for i = 1:nLong
    outputData{1,i} = longFields(1,i);
end
%put headers on the csv file for the modifier fields
for i = 1:nMod
    outputData{1,i+nLong} = modifierFields(1,i);
end

%iterate over the size of uLengths (the number of fields in the longest
%array
%keep track of which row in the 'flattened' data structure we are in
outR = 2;
for i = 1:size(dataStructure, 2)
    %pull out the length of one of the long fields
    firstLong = dataStructure(i).(longFields(1));
    nEntries = numel(firstLong);
    for j = 1:nEntries
        %write the value from the long (unique) arrays into the data
        %structure
        for k = 1:nLong
            arrayToInput = dataStructure(i).(longFields(k));
            %find entry j in the relevant long array
            outputData{outR,k} = arrayToInput(j);
        end
        %write all the modifier values that correspond to this value
        for k = 1:nMod
            outputData{outR,k+nLong} = dataStructure(i).(modifierFields(k));          
        end
        outR = outR + 1;
    end
end

%build the save path for the final version
fsep = filesep;
fullOutPath = strcat(outputPath, fsep, outputName ,'.csv');
%convert the cell array to a table for easier writing to a .csv file
tOut = cell2table(outputData);
%write the table tOut to the designated file path
writetable(tOut,fullOutPath,'writevariablenames',0);

end