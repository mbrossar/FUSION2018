function [trackerBis,myTracks] = razTrackerBis(ParamGlobal,ParamFilter,I,IdxImage,myTracks)
% reset trackerBis when number of points tracked is too small

fileImages = ParamGlobal.fileImages;
cameraParams = ParamFilter.cameraParams;
dirImage = ParamGlobal.dirImage;
idxImage = IdxImage-(0:1:20); %trakcing on previous images
idxReal = I-(0:10:10*length(idxImage)-10);
image = strcat(dirImage,int2str(fileImages(idxImage(1))),'.png');
image = undistortImage(imread(image),cameraParams);

trackerBis = vision.PointTracker('MaxBidirectionalError',1);
points = detectMinEigenFeatures(image);
points = selectUniform(points,200,size(image));
Location = points.Location;
initialize(trackerBis,Location,image);
trackerInit = clone(trackerBis);
for i = 1:length(Location)
    pT = pointTrack(1,Location(i,:));
    myTracks(i) = pT;
    myTracks(i).ViewIds = idxReal(1);
end
for i = 2:length(idxImage)
    image = strcat(dirImage,int2str(fileImages(idxImage(i))),'.png');
    image = undistortImage(imread(image),cameraParams);
    [points, validity] = trackerInit.step(image);
    for ii = 1:length(validity)
        if(validity(ii) == 1)
            myTracks(ii).ViewIds = [idxReal(i) myTracks(ii).ViewIds];
            myTracks(ii).Points = [points(ii,:); myTracks(ii).Points];
        end
    end
end
end

