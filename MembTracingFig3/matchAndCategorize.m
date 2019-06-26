%This function takes two strings that correspond to file names for FLIM and
%photon images as an input. It checks that the names match each other, and
%it returns the base name (without a file extension or a _color coded
%value/_photon modifier) to the calling function. 

%The function parses the file name strings to determine a coverslip and
%image number for each file. It assumes that the file names start with the
%format YYYY-MM-DD_XX-YY where XX is the coverslip ID and YY is the image
%ID. It also tolerates single digit image and coverslip IDs.

%Based on a metadata object that is passed to the function, it also
%categorizes each image. The metadata object should have two columns, with
%the first column to a unique coverslip ID and the second
%coverslip corresponding to the name of a category. Category names may (and
%should) be repeated.

%Last edited August 28, 2018 by Julia Lazzari-Dean

function [baseName, category, coverslipID, imageID] = matchAndCategorize(photonName,flimName,metadata)

%initialize the baseName to a nonsensical value for later error checking
baseName = 'error';

%determine the base file name for the photon images
baseP = char(photonName);
baseP = baseP(1:end-12);

%determine the base file name for the lifetime images
baseF = char(flimName);
baseF = baseF(1:end-22);

%if the base file names match, return this as the baseName
if(strcmp(baseF,baseP))
    baseName = baseF;
else
    error('Improperly paired photon and FLIM files.');
end

%start parsing the base file name to identify the coverslip and image IDs
%remove the date from the name to get to the coverslip ID
noDate = baseName(12:end);
%initialize coverslip ID and image IDs to a nonsensical value
coverslipID = -1;
imageID = -1;
imageIDExist = true(1);

%check if the coverslip ID is one or two digits by looking for the dash or
%underscore after it to specify the image number
if(noDate(3) == '-')
    %format is XX-XX
    coverslipID = str2double(noDate(1:2));
    updatedString = noDate(4:end);
elseif (noDate(2) == '-')
    %format is X-XX
    coverslipID = str2double(noDate(1));
    updatedString = noDate(3:end);
elseif(noDate(3) == '_')
    %format is XX_ (no imageID)
    coverslipID = str2double(noDate(1:2));
    imageIDExist = false;
elseif(noDate(2) == '_')
    %format is X_ (no imageID)
    imageIDExist = false;    
else
    error('Improperly formatted filename - cannot get coverslipID.')
end

%read further into the string to get the imageID
%check if the imageID is one or two digits by looking for the dash or
%underscore after it to specify the image number
if(updatedString(3) == '_' && imageIDExist)
    %two digit cimageID
    imageID = str2double(updatedString(1:2));
elseif (updatedString(2) == '_')
    %one digit imageID
    imageID = str2double(updatedString(1));
elseif (~imageIDExist)
    disp(['coverslip ID without imageID found for ' baseName]);
    imageID = 0;
else
    error('Improperly formatted filename - cannot get imageID.')
end

%match the category to a filename in the list
%set the category to a nonsensical value
category = 'error';
for i = 1:size(metadata,1)
    if(metadata{i,1} == coverslipID)
       category =  metadata{i,2};
    end
end
%return an error if the category was not found
if(strcmp(category,'error'))
    error('Could not categorize coverslip based on metadata.');
end


end