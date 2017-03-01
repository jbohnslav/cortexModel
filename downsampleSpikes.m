% spikes is a zero or one array 
% binsize is in milliseconds
% dt is in milliseconds
function [downsampledSpikes,downsampledT] = ...
    downsampleSpikes(spikes,binSize, dt)

binSize = binSize/dt;
NT = length(spikes);

binStarts = 1:binSize:NT;
binEnds = binSize:binSize:NT;
nBins = length(binEnds);
downsampledSpikes = nan(1,nBins);
for i=1:nBins
    downsampledSpikes(i) = sum(spikes(1,binStarts(i):binEnds(i)));
end
downsampledT = binStarts*dt;

end

