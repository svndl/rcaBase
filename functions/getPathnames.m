function [pathnames] = getPathnames(dataLocation,dataPrefix,dataSuffix)
% [pathnames] = getPathnames(dataLocation,dataPrefix,dataSuffix)
% 
% Prepares cell array of full paths to the text file locations to be read
% in by rcaSweep.m.
%
% If dataLocation is a cell array of different paths, it will loop through
% all and prepare one long concatenated list of pathnames. 
%
% If dataLocation is simply a string, the list of pathnames will be all
% matching files in dataLocation.


if ~iscell(dataLocation)
    folders = dir([dataLocation,'/',dataPrefix,'*',dataSuffix]);
    if isempty(folders)
        error('No files could be found for %s.\n',[dataLocation,'/',dataPrefix,'*',dataSuffix]);
    else
        numFolders = length(folders);
        pathnames = cell(1,numFolders);
        for k = 1:numFolders
            pathnames{k} = [dataLocation,'/',folders(k).name];
        end
    end
else
    numLocations = length(dataLocation);
    pathnames = {};
    for f = 1:numLocations
        folders = dir([dataLocation{f},'/',dataPrefix,'*',dataSuffix]);
        if isempty(folders)
            error('No files could be found for %s.\n',[dataLocation{f},'/',dataPrefix,'*',dataSuffix]);
        else
            numFolders = length(folders);
            prevLength = length(pathnames);
            for k = 1:numFolders
                pathnames{prevLength+k} = [dataLocation{f},'/',folders(k).name];
            end
        end
    end
end