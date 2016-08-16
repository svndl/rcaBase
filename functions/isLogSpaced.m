function [tIsLog] = isLogSpaced(binLevels)
% [tIsLog] = isLogSpaced(binLevels)

tStartBin = 1;
tNBins = length( binLevels );
tEndBin = tNBins;
tStart = binLevels( tStartBin );
tEnd = binLevels( tEndBin );

if( sum( abs( binLevels - logspace( log10(tStart), log10(tEnd), tNBins ) ) ) ...
        < sum( abs( binLevels - linspace( tStart, tEnd, tNBins ) ) ) )
    tIsLog = true;
else
    tIsLog = false;
end