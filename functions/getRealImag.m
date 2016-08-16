function [rcaDataReal,rcaDataImag] = getRealImag(rcaData)
% [rcaReal,rcaImag] = getRealImag(rcaData)
%   decompose rcaData for SSVEP paradigms into real and imaginary parts
%
%   07/20/2014. HEG. Note: will this work for non-sweep SSVEP rcaData?

% rcaDataReal and rcaDataImag are cell arrays of size
% nConditions-by-nSubjects, where each element of the cell array has
% dimensions coefficient-by-component-by-trial
%
% NB: the coefficients are arranged first across bins, and then
% concatenated across frequencies.  For example, if one has 10 bins and 2 frequencies, then the
% 20-element coefficient vector refers to
% [(bin1,f1),(bin2,f1),...,(bin10,f1),(bin1,f2),...(bin10,f2)]

rcaDataReal=cell(size(rcaData));
rcaDataImag=cell(size(rcaData));

% loop through elements of cell array and cut into real and imag halves
for i=1:size(rcaDataReal,1)
    for j=1:size(rcaDataReal,2)
   
        % check number of coefficients in this cell
        thisData=rcaData{i,j};
        nSamples=size(thisData,1);
        if rem(nSamples,2)~=0, error('odd number of samples is not consistent with ssvep data');  end
        
        % cut and store
        rcaDataReal{i,j}=rcaData{i,j}(1:nSamples/2,:,:);
        rcaDataImag{i,j}=rcaData{i,j}(nSamples/2+1:nSamples,:,:);
        
    end
end