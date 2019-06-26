%This function takes an MxN array inputData, where the columns N each
%represent unique experimental conditions. The function calculates the
%count, mean (meanTau - usually will be used with lifetime), median
%(medianTau), standard deviation (stdev), and standard error of the mean
%(sem) for each column. The function returns a struct outputStats, in which
%each column corresponds to one of the columns in the original array.

function [outputStats] = calculateStats(inputData)

%declare a structure to hold the output values
outputStats = struct('column',0,'count',0,'meanTau',0,'medianTau',0,'stdev',0,'sem',0);

%iterate across the input data and calculate the given stats for each
%column
for i = 1:size(inputData,2)
    outputStats(i).column = i;
    outputStats(i).count = sum(~isnan(inputData(:,i)));
    outputStats(i).meanTau = nanmean(inputData(:,i));
    outputStats(i).medianTau = nanmedian(inputData(:,i));
    outputStats(i).stdev = nanstd(inputData(:,i));
    outputStats(i).sem = outputStats(i).stdev/sqrt(outputStats(i).count);
end

end
