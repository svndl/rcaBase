function cmplxData= getCmplx(realPart, imagPart)
%  cmplxData= getCmplx(realPart, imagPart)
%   join real and imaginary parts of data (e.g. for estimaing rca or ssd) for SSVEP paradigms into complex representation
%
%   07/20/2014. HEG. Note: will this work for non-sweep SSVEP rcaData?

% realPart, imagPart and cmplxData are cell arrays of size
% nConditions-by-nSubjects, where each element of the cell array has
% dimensions coefficient-by-component-by-trial
%
% NB: the coefficients are arranged first across bins, and then
% concatenated across frequencies.  For example, if one has 10 bins and 2 frequencies, then the
% 20-element coefficient vector refers to
% [(bin1,f1),(bin2,f1),...,(bin10,f1),(bin1,f2),...(bin10,f2)]
% this follows exactly the nomenclature of getRealImag.m

cmplxData =cell(size(realPart));

% loop through elements of cell array and cut into real and imag halves
for i=1:size(cmplxData ,1)
    for j=1:size(cmplxData ,2)
   
        % check number of coefficients in this cell
        thisRealPart = realPart{i,j};
        thisImagPart = imagPart{i,j};
        
        nSamples=size(thisRealPart,1);
        if size(thisRealPart)~=size(thisImagPart), error('real and imag part of data do not have the same size'); end
        
        % cut and store
        cmplxData{i,j}  = complex(thisRealPart, thisImagPart) ; 
    end
end