
function [ ind] = feGetRegionVox(fe, regions)

ind = [];
for i=1:length(regions)
    ind = [ind; find(fe.life.M.segmentation == regions(i))];
end

end

